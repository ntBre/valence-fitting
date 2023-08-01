GENERATE := 01_generate-forcefield
CURATE := 02_curate-data
MSM := 03_generate-initial-ff
FIT := 04_fit-forcefield

INITIAL_FF := $(GENERATE)/output/initial-force-field-openff-2.1.0.offxml
OPT_SET := $(addprefix $(CURATE)/output/,pavan-opt-training-set.json pavan-opt-smirks.json)
TD_SET := $(addprefix $(CURATE)/output/,pavan-td-training-set.json pavan-td-smirks.json)
COMBINED := $(addprefix $(CURATE)/output/,combined-opt.json combined-opt-smirks.json combined-td.json combined-td-smirks.json)
MSM_FF := $(MSM)/output/initial-force-field-msm.offxml

.PHONY: step1 step2 step3 step4 sage td

step1: $(INITIAL_FF)
step2: $(COMBINED)
step3: $(MSM_FF)
step4: $(FIT)/ready
sage: fit-sage/ready

td: $(TD_SET)
opt: $(OPT_SET)

$(INITIAL_FF): $(GENERATE)/generate-forcefield.py
	echo $$OE_LICENSE
	cd $(GENERATE) ; \
	python generate-forcefield.py --force-field-name openff-2.1.0.offxml \
                                      --output           output/initial-force-field-openff-2.1.0.offxml \
                                      2>&1 | tee generate.log

$(OPT_SET) &: $(INITIAL_FF) $(CURATE)/curate_dataset.py
	cd $(CURATE) ; \
	python curate_dataset.py download-opt                                          \
    --core-opt-dataset          "OpenFF multiplicity correction optimization set v1.0" \
    --initial-forcefield        ../$(INITIAL_FF)                                       \
    --max-opt-conformers        12                                                     \
    --n-processes               8                                                      \
    --min-record-coverage       1                                                      \
    --opt-records-to-remove     opt_records_to_remove.dat                              \
    --output                    "output/pavan-opt-training-set.json"                   \
    --output-parameter-smirks   "output/pavan-opt-smirks.json"                         \
    --verbose

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
    --verbose ; \
	date > ready

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
