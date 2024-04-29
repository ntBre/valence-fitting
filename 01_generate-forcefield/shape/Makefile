ff := ../output/initial-force-field-openff-2.1.0.offxml

labels.json: label_params.py
	python $< -i datasets/qm/output -f $(ff) -o $@

ff-plots:
	python plot_ff_torsions.py -f $(ff) -o figs

td-plots:
	python plot_tds.py \
-f $(ff) \
-q datasets/qm/output/torsiondrive \
-l labels.json \
-p t60g
