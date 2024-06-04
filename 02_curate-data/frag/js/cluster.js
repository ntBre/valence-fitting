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

function drawMolecule(moldata) {
	let dialog = document.getElementById("edit-molecule-modal");

	// new frame to put into the dialog
	let frame = document.createElement("div");
	frame.setAttribute("id", "edit-molecule-modal-content");

	let canvas = document.createElement("canvas");
	canvas.setAttribute("width", moldata.canvas_size);
	canvas.setAttribute("height", moldata.canvas_size);
	let ctx = canvas.getContext("2d");
	ctx.fillStyle = "black";
	for (let i = 0; i < moldata.atoms.length; i++) {
		let x = moldata.coords[i][0];
		let y = moldata.coords[i][1];
		let atomic_sym = moldata.atoms[i];
		if (atomic_sym != "C") {
			ctx.fillText(atomic_sym, x, y);
		}
	}

	for ([a1, a2, order] of moldata.bonds) {
		let hl_bond = moldata.hl_atoms.includes(a1) && moldata.hl_atoms.includes(a2);
		let [x1, y1] = moldata.coords[a1];
		let [x2, y2] = moldata.coords[a2];
		if (hl_bond) {
			ctx.strokeStyle = "orange";
		}
		if (order == 1) {
			ctx.beginPath();
			ctx.moveTo(x1, y1);
			ctx.lineTo(x2, y2);
			ctx.stroke();
		} else if (order === 2) {
			// compute perpendicular unit vector to move along for double bonds
			let dx = x2 - x1;
			let dy = y2 - y1;
			let m = Math.sqrt(dx * dx + dy * dy);
			let [ux, uy] = [-dy / m, dx / m];

			let f = 2;

			ctx.beginPath();
			ctx.moveTo(x1 + f * ux, y1 + f * uy);
			ctx.lineTo(x2 + f * ux, y2 + f * uy);
			ctx.stroke();

			ctx.beginPath();
			ctx.moveTo(x1 - f * ux, y1 - f * uy);
			ctx.lineTo(x2 - f * ux, y2 - f * uy);
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
