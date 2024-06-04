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
		console.log(moldata);
		drawMolecule(moldata);
	} else {
		console.log("no smiles found");
	}
}

function drawMolecule(moldata) {
	console.log("in drawMolecule");
	let dialog = document.getElementById("edit-molecule-modal");

	// new frame to put into the dialog
	let frame = document.createElement("div");
	frame.setAttribute("id", "edit-molecule-modal-content");

	let canvas = document.createElement("canvas");
	canvas.setAttribute("width", "400");
	canvas.setAttribute("height", "400");
	let ctx = canvas.getContext("2d");
	ctx.fillStyle = "black";
	for (let i = 0; i < moldata.atoms.length; i++) {
		let x = moldata.coords[i][0];
		let y = moldata.coords[i][1];
		ctx.arc(x, y, 3, 0, 2*Math.PI);
		ctx.fill();
		ctx.closePath();
	}

	for (b of moldata.bonds) {
		let c1 = moldata.coords[b[0]];
		let c2 = moldata.coords[b[1]];
		console.log(b, c1, c2);
		let x1 = c1[0];
		let y1 = c1[1];
		let x2 = c2[0];
		let y2 = c2[1];
		ctx.beginPath();
		ctx.moveTo(x1, y1);
		ctx.lineTo(x2, y2);
		ctx.stroke();
		console.log(b, moldata.atoms[b[0]], moldata.atoms[b[1]]);
	}

	frame.appendChild(canvas);

	// do replace
	let toReplace = document.getElementById("edit-molecule-modal-content");
	dialog.replaceChild(frame, toReplace);
	dialog.showModal();
}
