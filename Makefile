GENERATE := 01_generate-forcefield
CURATE := 02_curate-data
MSM := 03_generate-initial-ff
FIT := 04_fit-forcefield

MSM_FF := $(MSM)/output/initial-force-field-msm.offxml

.PHONY: step1 step2 step3 step4 sage td opt

step3: $(MSM_FF)
sage: fit-sage/ready

# step 1 - generate the initial force field from a combination of the
# differences between Sage 2.0 and pavan's force field and Sage 2.1 (port the
# changes Pavan made to Sage 2.1). Inputs are the script and Pavan's FF, output
# is the initial force field

INITIAL_FF := $(GENERATE)/output/initial-force-field-openff-2.1.0.offxml
STEP1_DEPS := $(GENERATE)/generate-forcefield.py $(GENERATE)/force-field.offxml
$(INITIAL_FF): $(STEP1_DEPS)
	cd $(GENERATE); python generate-forcefield.py --output ../$@

step1: $(INITIAL_FF) # alias for step1

# step 2 - curate data sets for fitting the force field. In this case, we want
# to combine the original data sets used for training Sage 2.1 with some
# additional records selected by Pavan specifically for the additional torsions
# added to the force field.

# step 2a - download, filter, and curate the torsion multiplicity optimization
# data set. Inputs are
# 1) The initial force field from step 1
# 2) curate_dataset.py
# 3) filtered-core-opt.json (pre-filtered opt dataset)
# Outputs are
# 1) pavan-opt-training-set.json (the final data set)
# 2) pavan-opt-smirks.json (final smirks)

# step 2a.i - download the data set with qcportal. This only depends on the
# script itself and outputs the dataset
$(CURATE)/datasets/core-opt.json: $(CURATE)/download_opt.py
	cd $(CURATE); python download_opt.py

# step 2a.ii - filter the opt data set with fast-filter. This depends on the
# filtering code in filters.py, the script passed to fast-filter (filter-opt.py)
# and the downloaded data set from 2a.i. The output is the filtered dataset
DEPS := $(addprefix $(CURATE)/,datasets/core-opt.json filters.py filter-opt.py)
$(CURATE)/datasets/filtered-core-opt.json: $(DEPS)
	cd $(CURATE); \
	fast-filter datasets/core-opt.json -p filter-opt.py -o ../$@ -t 12

# step 2a.iii - perform the last steps of curation - just filtering duplicate
# records I guess, and then select_parameters for the smirks output. This simply
# reads in the filtered data set, so it depends on that and curate_dataset.py.
# It also uses the initial FF for select_parameters. Outputs are the training
# set and smirks
OPT_SET := $(addprefix $(CURATE)/output/,pavan-opt-training-set.json	\
pavan-opt-smirks.json)

DEPS := $(INITIAL_FF) $(CURATE)/curate_dataset.py	\
$(CURATE)/datasets/filtered-core-opt.json

$(OPT_SET) &: $(DEPS)
	cd $(CURATE) ; \
	python curate_dataset.py download-opt                                \
	    --initial-forcefield        ../$(INITIAL_FF)                     \
	    --min-record-coverage       1                                    \
	    --output                    "output/pavan-opt-training-set.json" \
	    --output-parameter-smirks   "output/pavan-opt-smirks.json"       \
	    --verbose

opt: $(OPT_SET)

# step 2b - the same as 2a but for the torsion data sets

# step 2b.i - download the data set with qcportal
$(CURATE)/datasets/core-td.json: $(CURATE)/download_td.py
	cd $(CURATE); python download_td.py

# step 2b.ii - filter the td data set with fast-filter
DEPS := $(addprefix $(CURATE)/,datasets/core-td.json filters.py filter-td.py)
$(CURATE)/datasets/filtered-core-td.json: $(DEPS)
	cd $(CURATE); \
	fast-filter datasets/core-td.json -p filter-td.py -o ../$@ -t 12

# step 2b.iii - last steps of curation
TD_SET := $(addprefix $(CURATE)/output/,pavan-td-training-set.json	\
	    pavan-td-smirks.json)

DEPS := $(INITIAL_FF) $(CURATE)/curate_dataset.py	\
	$(CURATE)/datasets/filtered-core-td.json

$(TD_SET) &: $(DEPS)
	cd $(CURATE) ; \
	python curate_dataset.py download-td                                \
	    --initial-forcefield        ../$(INITIAL_FF)                    \
	    --output                    "output/pavan-td-training-set.json" \
	    --output-parameter-smirks   "output/pavan-td-smirks.json"       \
	    --n-processes               8                                   \
	    --min-record-coverage       1                                   \
	    --verbose

td: $(TD_SET)

# step 2c curate the sage datasets using the same filters as above

# step 2c.i filter the sage opt dataset
DEPS := $(addprefix $(CURATE)/,sage/opt-set-for-fitting-2.1.0.json	\
	filter-opt.py filters.py)

$(CURATE)/datasets/filtered-sage-opt.json: $(DEPS)
	cd $(CURATE) ; \
	fast-filter sage/opt-set-for-fitting-2.1.0.json \
		-p filter-opt.py -o ../$@ -t 12 -b 32

# step 2c.ii filter the sage td dataset
DEPS := $(addprefix $(CURATE)/,sage/td-set-for-fitting-2.1.0.json filter-td.py	\
	filters.py)

$(CURATE)/datasets/filtered-sage-td.json: $(DEPS)
	cd $(CURATE) ; \
	fast-filter sage/td-set-for-fitting-2.1.0.json \
		-p filter-td.py -o ../$@ -t 12 -b 32

DEPS := $(OPT_SET) $(TD_SET) $(CURATE)/combine.py	\
	$(CURATE)/datasets/filtered-sage-opt.json	\
	$(CURATE)/datasets/filtered-sage-td.json

# step 2d combine the two pairs of filtered data sets together
COMBINED := $(addprefix $(CURATE)/output/,combined-opt.json			\
	    combined-opt-smirks.json combined-td.json combined-td-smirks.json)

$(COMBINED) &: $(DEPS)
	cd $(CURATE); python combine.py

step2: $(COMBINED)



# step 3 - update the initial force field parameters with the modified seminario
# method

$(MSM_FF): $(INITIAL_FF) $(CURATE)/output/combined-opt.json $(MSM)/create-msm-ff.py
	cd $(MSM) ; \
	python create-msm-ff.py                                                \
	--initial-force-field       ../$(INITIAL_FF)                           \
	--optimization-dataset      ../02_curate-data/output/combined-opt.json \
	--working-directory         working-directory                          \
	--output                    ../$@

# step 4 - generate ForceBalance inputs

DEPS := $(FIT)/smiles-to-exclude.dat $(FIT)/smarts-to-exclude.dat
$(FIT)/ready: $(COMBINED) $(MSM_FF) $(FIT)/create-fb-inputs.py $(DEPS)
	rm -r $(FIT)/fb-fit/targets
	mkdir -p $(FIT)/fb-fit/targets
	cd $(FIT) ; \
	python create-fb-inputs.py                                                      \
	--tag                       "fb-fit"                                            \
	--optimization-dataset      "../02_curate-data/output/combined-opt.json"        \
	--torsion-dataset           "../02_curate-data/output/combined-td.json"         \
	--forcefield                ../$(MSM_FF)                                        \
	--valence-to-optimize       "../02_curate-data/output/combined-opt-smirks.json" \
	--torsions-to-optimize      "../02_curate-data/output/combined-td-smirks.json"  \
	--smiles-to-exclude         "smiles-to-exclude.dat"                             \
	--smarts-to-exclude         "smarts-to-exclude.dat"                             \
	--max-iterations            1                                                   \
	--port                      55387                                               \
	--output-directory          "output"                                            \
	--verbose
	date > $@

step4: $(FIT)/ready

# step 5 - pack up generated files for running on HPC3

# step 5a - pack up the targets directory
$(FIT)/fb-fit/targets.tar.gz: $(wildcard $(FIT)/fb-fit/targets/*/*)
	cd $(FIT)/fb-fit ; tar cfz targets.tar.gz targets

# step 5b - zip everything up
TORS_DEPS := $(addprefix $(FIT)/,$(addprefix					\
				        fb-fit/,forcefield/force-field.offxml	\
				        optimize.in targets.tar.gz))

tors.tar.gz: $(FIT)/ready $(TORS_DEPS)
	tar cfz $@ $(TORS_DEPS)
