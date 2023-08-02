GENERATE := 01_generate-forcefield
CURATE := 02_curate-data
MSM := 03_generate-initial-ff
FIT := 04_fit-forcefield

TD_SET := $(addprefix $(CURATE)/output/,pavan-td-training-set.json pavan-td-smirks.json)
COMBINED := $(addprefix $(CURATE)/output/,combined-opt.json combined-opt-smirks.json combined-td.json combined-td-smirks.json)
MSM_FF := $(MSM)/output/initial-force-field-msm.offxml

.PHONY: step1 step2 step3 step4 sage td opt

step2: $(COMBINED)
step3: $(MSM_FF)
step4: $(FIT)/ready
sage: fit-sage/ready

td: $(TD_SET)

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

# step 2a - download the torsion multiplicity optimization data set. Inputs are
# 1) The initial force field from step 1
# 2) The script itself
# Primary outputs are
# 1) pavan-opt-training-set.json (the final data set)
# 2) pavan-opt-smirks.json
# Intermediate outputs are
# 1) datasets/core-opt.json (downloaded core dataset)
# 2) datasets/filtered-core-opt.json (filtered core dataset)

# these caches are actually potential inputs. In fact, the second cache needs to
# be invalidated if there are any changes to the filter

$(CURATE)/datasets/core-opt.json: $(CURATE)/download_opt.py
	cd $(CURATE); python download_opt.py

DEPS := $(addprefix $(CURATE)/,datasets/core-opt.json filters.py filter-opt.py)
$(CURATE)/datasets/filtered-core-opt.json: $(DEPS)
	cd $(CURATE); \
	fast-filter datasets/core-opt.json -p filter-opt.py -o ../$@ -t 12

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

$(TD_SET) &: $(INITIAL_FF) $(CURATE)/curate_dataset.py
	cd $(CURATE) ; \
	python curate_dataset.py download-td                                             \
    --core-td-dataset           "OpenFF multiplicity correction torsion drive data v1.1" \
    --initial-forcefield        ../$(INITIAL_FF)                                         \
    --n-processes               8                                                        \
    --min-record-coverage       1                                                        \
    --td-records-to-remove      td_records_to_remove.dat                                 \
    --output                    "output/pavan-td-training-set.json"                      \
    --output-parameter-smirks   "output/pavan-td-smirks.json"                            \
    --verbose

# run fast-filter to generate the cache read by combine.py
SAGE_DATA_SETS := ../../clone/sage-2.1.0/inputs-and-outputs/data-sets
$(CURATE)/datasets/filtered-sage-opt.json: $(SAGE_DATA_SETS)/opt-set-for-fitting-2.1.0.json $(CURATE)/filter-opt.py
	cd $(CURATE) ; \
	fast-filter ../$(SAGE_DATA_SETS)/opt-set-for-fitting-2.1.0.json -p filter-opt.py -o ../$@ -t 12 -b 32

$(COMBINED) &: $(OPT_SET) $(TD_SET) $(CURATE)/combine.py $(CURATE)/datasets/filtered-sage-opt.json
	cd $(CURATE) ; \
	python combine.py

$(MSM_FF): $(INITIAL_FF) $(CURATE)/output/combined-opt.json $(MSM)/create-msm-ff.py
	cd $(MSM) ; \
	python create-msm-ff.py                                                                            \
    --initial-force-field       "../01_generate-forcefield/output/initial-force-field-openff-2.1.0.offxml" \
    --optimization-dataset      "../02_curate-data/output/combined-opt.json"                               \
    --working-directory         "working-directory"                                                        \
    --output                    "output/initial-force-field-msm.offxml"

$(FIT)/ready: $(COMBINED) $(MSM_FF) $(FIT)/create-fb-inputs.py
	cd $(FIT) ; \
	python create-fb-inputs.py                                                                \
	    --tag                       "fb-fit"                                                          \
	    --optimization-dataset      "../02_curate-data/output/combined-opt.json"                      \
	    --torsion-dataset           "../02_curate-data/output/combined-td.json"                       \
	    --forcefield                "../03_generate-initial-ff/output/initial-force-field-msm.offxml" \
	    --valence-to-optimize       "../02_curate-data/output/combined-opt-smirks.json"               \
	    --torsions-to-optimize      "../02_curate-data/output/combined-td-smirks.json"                \
	    --smiles-to-exclude         "smiles-to-exclude.dat"                                           \
	    --smarts-to-exclude         "smarts-to-exclude.dat"                                           \
	    --max-iterations            1                                                                 \
	    --port                      55387                                                             \
	    --output-directory          "output"                                                          \
	    --verbose
	date > $@

fit-sage/ready: $(COMBINED) sage-2.1.0-msm.offxml
	cd fit-sage ; \
	python ../$(FIT)/create-fb-inputs.py                                                      \
    --tag                       "fb-fit"                                                          \
    --optimization-dataset      "../02_curate-data/output/combined-opt.json"                      \
    --torsion-dataset           "../02_curate-data/output/combined-td.json"                       \
    --valence-to-optimize       "../02_curate-data/output/combined-opt-smirks.json"               \
    --torsions-to-optimize      "../02_curate-data/output/combined-td-smirks.json"                \
    --forcefield                "../sage-2.1.0-msm.offxml"                                        \
    --smiles-to-exclude         ../$(FIT)/smiles-to-exclude.dat                                   \
    --smarts-to-exclude         ../$(FIT)/smarts-to-exclude.dat                                   \
    --max-iterations            1                                                                 \
    --port                      55387                                                             \
    --output-directory          "output"                                                          \
    --verbose
	date > $@

DEPS :=$(addprefix fit-sage/, $(addprefix fb-fit/,forcefield/force-field.offxml	\
       optimize.in targets.tar.gz) parameters-to-optimize scripts)

fit-sage/fb-fit/targets.tar.gz: $(wildcard fit-sage/fb-fit/targets/*)
	cd fit-sage/fb-fit ; \
	tar cvfz targets.tar.gz targets

sage.tar.gz: $(DEPS)
	tar cvfz $@ $^

$(FIT)/fb-fit/targets.tar.gz: $(wildcard $(FIT)/fb-fit/targets/*/*)
	cd $(FIT)/fb-fit ; \
	tar cvfz targets.tar.gz targets

TORS_DEPS := $(addprefix $(FIT)/,$(addprefix fb-fit/,forcefield/force-field.offxml	\
				        optimize.in targets.tar.gz))

tors.tar.gz: $(FIT)/ready $(FIT)/fb-fit/targets.tar.gz
	tar cvfz $@ $(TORS_DEPS)
