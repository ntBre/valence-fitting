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
	for (let i = 0; i < mol.atoms.length; i++) {
		let x = mol.coords[i][0];
		let y = mol.coords[i][1];
		let [atomic_sym, charge] = mol.atoms[i];
		if (charge != 0 || atomic_sym != "C") {
			if (charge == 0) {
				centerText(ctx, atomic_sym, x, y, font_size);
			} else if (charge == -1) {
				centerText(ctx, atomic_sym + "-", x, y, font_size);
			} else if (charge == 1) {
				centerText(ctx, atomic_sym + "+", x, y, font_size);
			} else {
				console.log("warning: unrecognized atomic charge: ", charge);
			}
		}
	}

	for ([a1, a2, order] of mol.bonds) {
		let hl_bond = mol.hl_atoms.includes(a1) && mol.hl_atoms.includes(a2);
		let [x1, y1] = mol.coords[a1];
		let [x2, y2] = mol.coords[a2];
		let dx = x2 - x1;
		let dy = y2 - y1;
		let m = Math.sqrt(dx * dx + dy * dy);
		let [ux, uy] = [dx / m, dy / m];

		let [as1, ch1] = mol.atoms[a1];
		let [as2, ch2] = mol.atoms[a2];
		let pad = 8;
		if (ch1 != 0 || as1 != "C") {
			x1 += pad * ux;
			y1 += pad * uy;
		}
		if (ch2 != 0 || as2 != "C") {
			x2 -= pad * ux;
			y2 -= pad * uy;
		}
		// perpendicular unit vector to move along for double and triple bonds
		let [px, py] = [-uy, ux];
		if (hl_bond) {
			ctx.strokeStyle = "orange";
		}
		if (order == 1) {
			ctx.beginPath();
			ctx.moveTo(x1, y1);
			ctx.lineTo(x2, y2);
			ctx.stroke();
		} else if (order === 2) {
			let f = 2;

			ctx.beginPath();
			ctx.moveTo(x1 + f * px, y1 + f * py);
			ctx.lineTo(x2 + f * px, y2 + f * py);
			ctx.stroke();

			ctx.beginPath();
			ctx.moveTo(x1 - f * px, y1 - f * py);
			ctx.lineTo(x2 - f * px, y2 - f * py);
			ctx.stroke();
		} else {
			console.log("warning: unknown bond order ", order);
			ctx.strokeStyle = "red";
			ctx.beginPath();
			ctx.moveTo(x1, y1);
			ctx.lineTo(x2, y2);
			ctx.stroke();
		}
		ctx.strokeStyle = "black";
	}

	frame.appendChild(canvas);

	// do replace
	let toReplace = document.getElementById("edit-molecule-modal-content");
	dialog.replaceChild(frame, toReplace);
	dialog.showModal();
}
