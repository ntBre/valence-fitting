async function addToDataset() {
	let node = document.getElementById("modal-box-content");
	let smiles = node.getAttribute("smiles");
	let pid = node.getAttribute("pid");
	if (smiles) {
		let response = await fetch("/add-molecule", {
			method: "POST",
			headers: {
				'Accept': 'application/json',
				'Content-Type': 'application/json'
			},
			body: JSON.stringify({
				"smiles": smiles,
				"pid": pid
			}),
		});
		if (!response.ok) {
			console.log(`error handling request: ${response}`);
		}
	} else {
		console.log("no smiles found");
	}
}

async function editMolecule() {
	let node = document.getElementById("modal-box-content");
	let smiles = node.getAttribute("smiles");
	let pid = node.getAttribute("pid");
	if (smiles) {
		let response = await fetch("/edit-molecule", {
			method: "POST",
			headers: {
				'Accept': 'application/json',
				'Content-Type': 'application/json'
			},
			body: JSON.stringify({
				"smiles": smiles,
				"pid": pid
			}),
		});
		if (!response.ok) {
			console.log(`error handling request: ${response}`);
			return;
		}
		let moldata = await response.json();
		drawMolecule(moldata);
	} else {
		console.log("no smiles found");
	}
}

function centerText(ctx, text, x, y, font_size) {
	let w = ctx.measureText(text).width;
	ctx.fillText(text, x - w / 2, y + font_size / 2);
}

class Atom {
	constructor(x, y, symbol, charge) {
		this.x = x;
		this.y = y;
		if (charge == 0) {
			this.symbol = symbol;
		} else if (charge == -1) {
			this.symbol = symbol + "-"
		} else if (charge == 1) {
			this.symbol = symbol + "+";
		} else {
			console.log("warning: unrecognized atomic charge: ", charge);
		}
		this.is_highlighted = false;
		this.is_selected = false;
	}

	draw(ctx, font_size) {
		if (this.symbol !== "C") {
			let w = ctx.measureText(this.symbol).width;
			ctx.fillText(this.symbol, this.x - w / 2, this.y + font_size / 2);
		}
		if (this.is_selected) {
			ctx.strokeStyle = "blue";
			ctx.beginPath();
			ctx.arc(this.x, this.y, 1.0 * font_size, 0, 2 * Math.PI);
			ctx.stroke();
			ctx.strokeStyle = "black";
		} else if (this.is_highlighted) {
			ctx.beginPath();
			ctx.arc(this.x, this.y, 1.0 * font_size, 0, 2 * Math.PI);
			ctx.stroke();
		}
	}

	contains(x, y, r) {
		let dx = x - this.x;
		let dy = y - this.y;
		return dx * dx + dy * dy <= r * r;
	}
}

class Bond {
	constructor(x1, y1, x2, y2, order, hl_bond) {
		let dx = x2 - x1;
		let dy = y2 - y1;
		let m = Math.sqrt(dx * dx + dy * dy);
		let [ux, uy] = [dx / m, dy / m];
		let [px, py] = [-uy, ux];
		this.x1 = x1;
		this.x2 = x2;
		this.y1 = y1;
		this.y2 = y2;
		this.ux = ux;
		this.uy = uy;
		this.px = px;
		this.py = py;
		this.order = order;
		this.is_highlighted = false;
		this.is_selected = false;
		this.hl_bond = hl_bond;
		let [cx, cy] = midpoint(x1, y1, x2, y2);
		this.midpoint = [cx + x1, cy + y1];
	}

	contains(x, y, r) {
		let [cx, cy] = this.midpoint;
		let dx = x - cx;
		let dy = y - cy;
		return dx * dx + dy * dy <= r * r;
	}

	draw(ctx, font_size) {
		const f = 2;
		if (this.hl_bond) {
			ctx.strokeStyle = "orange";
		}
		if (this.order == 1) {
			ctx.beginPath();
			ctx.moveTo(this.x1, this.y1);
			ctx.lineTo(this.x2, this.y2);
			ctx.stroke();
		} else if (this.order === 2) {
			ctx.beginPath();
			ctx.moveTo(this.x1 + f * this.px, this.y1 + f * this.py);
			ctx.lineTo(this.x2 + f * this.px, this.y2 + f * this.py);
			ctx.stroke();

			ctx.beginPath();
			ctx.moveTo(this.x1 - f * this.px, this.y1 - f * this.py);
			ctx.lineTo(this.x2 - f * this.px, this.y2 - f * this.py);
			ctx.stroke();
		} else {
			console.log("warning: unknown bond order ", this.order);
			ctx.strokeStyle = "red";
			ctx.beginPath();
			ctx.moveTo(this.x1, this.y1);
			ctx.lineTo(this.x2, this.y2);
			ctx.stroke();
		}
		ctx.strokeStyle = "black";
		if (this.is_selected) {
			let [cx, cy] = this.midpoint;
			ctx.strokeStyle = "blue";
			ctx.beginPath();
			ctx.arc(cx, cy, 0.5 * font_size, 0, 2 * Math.PI);
			ctx.stroke();
			ctx.strokeStyle = "black";
		} else if (this.is_highlighted) {
			let [cx, cy] = this.midpoint;
			ctx.beginPath();
			ctx.arc(cx, cy, 0.5 * font_size, 0, 2 * Math.PI);
			ctx.stroke();
		}
	}
}

function midpoint(x1, y1, x2, y2) {
	return [(x2 - x1) / 2, (y2 - y1) / 2];
}

class Scene {
	constructor(mol) {
		this.atoms = Array();
		for (let i = 0; i < mol.atoms.length; i++) {
			let x = mol.coords[i][0];
			let y = mol.coords[i][1];
			let [atomic_sym, charge] = mol.atoms[i];
			this.atoms.push(new Atom(x, y, atomic_sym, charge));
		}

		this.bonds = Array();
		let pad = 8; // padding for ends of bonds
		for (let [a1, a2, order] of mol.bonds) {
			let hl_bond = mol.hl_atoms.includes(a1) && mol.hl_atoms.includes(a2);
			let [x1, y1] = mol.coords[a1];
			let [x2, y2] = mol.coords[a2];
			let dx = x2 - x1;
			let dy = y2 - y1;
			let m = Math.sqrt(dx * dx + dy * dy);
			let [ux, uy] = [dx / m, dy / m];

			let [as1, ch1] = mol.atoms[a1];
			let [as2, ch2] = mol.atoms[a2];
			if (ch1 != 0 || as1 != "C") {
				x1 += pad * ux;
				y1 += pad * uy;
			}
			if (ch2 != 0 || as2 != "C") {
				x2 -= pad * ux;
				y2 -= pad * uy;
			}
			this.bonds.push(new Bond(x1, y1, x2, y2, order, hl_bond));
		}
	}

	draw(ctx, font_size) {
		for (let atom of this.atoms) {
			atom.draw(ctx, font_size);
		}
		for (let bond of this.bonds) {
			bond.draw(ctx, font_size);
		}
	}
}


function drawMolecule(mol) {
	let dialog = document.getElementById("edit-molecule-modal");

	// new frame to put into the dialog
	let frame = document.createElement("div");
	frame.setAttribute("id", "edit-molecule-modal-content");

	let canvas = document.createElement("canvas");
	canvas.setAttribute("width", mol.canvas_size);
	canvas.setAttribute("height", mol.canvas_size);
	let ctx = canvas.getContext("2d");
	let font_size = 20;
	ctx.font = `normal ${font_size}px sans-serif`;
	ctx.fillStyle = "black";
	ctx.lineWidth = 2;

	let scene = new Scene(mol);

	scene.draw(ctx, font_size);

	canvas.addEventListener("mousemove", (event) => {
		ctx.clearRect(0, 0, canvas.width, canvas.height);
		for (let atom of scene.atoms) {
			atom.is_highlighted = atom.contains(event.offsetX, event.offsetY, font_size);
		}
		for (let bond of scene.bonds) {
			bond.is_highlighted = bond.contains(event.offsetX, event.offsetY, font_size);
		}
		scene.draw(ctx, font_size);
	});

	canvas.addEventListener("click", (event) => {
		for (let atom of scene.atoms) {
			if (atom.contains(event.offsetX, event.offsetY, font_size)) {
				atom.is_selected = !atom.is_selected;
			}
		}
		for (let bond of scene.bonds) {
			if (bond.contains(event.offsetX, event.offsetY, font_size)) {
				bond.is_selected = !bond.is_selected;
			}
		}
	});

	frame.appendChild(canvas);

	// do replace
	let toReplace = document.getElementById("edit-molecule-modal-content");
	dialog.replaceChild(frame, toReplace);
	dialog.showModal();
}
