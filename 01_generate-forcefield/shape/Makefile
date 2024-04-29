ff := ../output/initial-force-field-openff-2.1.0.offxml

labels.json: label_params.py
	python $< -i datasets/qm/output -f $(ff) -o $@

plots:
	python plot_ff_torsions.py -f $(ff) -o figs
