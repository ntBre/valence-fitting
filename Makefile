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
STEP1_DEPS := $(GENERATE)/generate-forcefield.py $(GENERATE)/to_add.dat
$(INITIAL_FF): $(STEP1_DEPS)
	cd $(GENERATE); python generate-forcefield.py --output ../$@

step1: $(INITIAL_FF) # alias for step1

# step 2 - curate data sets for fitting the force field. In this case, we want
# to combine the original data sets used for training Sage 2.1 with some
# additional records selected by Pavan specifically for the additional torsions
# added to the force field.

## step 2a download the core-opt and core-td datasets
$(CURATE)/datasets/core-opt.json: $(CURATE)/download_opt.py
	cd $(CURATE); python download_opt.py

$(CURATE)/datasets/core-td.json: $(CURATE)/download_td.py
	cd $(CURATE); python download_td.py

$(CURATE)/datasets/supp-td.json: $(CURATE)/download_td.py
	cd $(CURATE); python download_td.py -o ../$@ -d "OpenFF Torsion Coverage Supplement v1.0"

$(CURATE)/datasets/supp2-td.json: $(CURATE)/download_td.py
	cd $(CURATE); python download_td.py -o ../$@ -d "OpenFF Torsion Multiplicity Torsion Drive Coverage Supplement v1.0"

$(CURATE)/datasets/supp-opt.json: $(CURATE)/download_opt.py
	cd $(CURATE); python $< -o ../$@ -d "OpenFF Torsion Multiplicity Optimization Training Coverage Supplement v1.0"

## step 2b filter sage data sets for charge issues
$(CURATE)/sage/filtered-opt.json: $(CURATE)/sage/opt.json $(CURATE)/charge-filter.py
	python $(CURATE)/charge-filter.py --input $< --output $@
$(CURATE)/sage/filtered-td.json: $(CURATE)/sage/td.json $(CURATE)/charge-filter.py
	python $(CURATE)/charge-filter.py --input $< --output $@

## step 2c fully filter the new datasets
DEPS := $(CURATE)/datasets/core-opt.json $(CURATE)/filters.py	\
	$(CURATE)/filter-opt.py

$(CURATE)/datasets/filtered-opt.json: $(DEPS)
	cd $(CURATE) ; python filter-opt.py --input ../$< --output ../$@

DEPS := $(CURATE)/datasets/core-td.json $(CURATE)/filters.py	\
	$(CURATE)/filter-td.py

$(CURATE)/datasets/filtered-td.json: $(DEPS)
	cd $(CURATE) ; python filter-td.py --input ../$< --output ../$@

$(CURATE)/datasets/filtered-%-td.json: $(CURATE)/datasets/%-td.json $(CURATE)/filters.py $(CURATE)/filter_supplement.py
	cd $(CURATE) ; python filter_supplement.py --input ../$< --output ../$@

$(CURATE)/datasets/filtered-%-opt.json: $(CURATE)/datasets/%-opt.json $(CURATE)/filters.py $(CURATE)/filter_opt_supplement.py
	cd $(CURATE) ; python filter_opt_supplement.py --input ../$< --output ../$@

## step 2d combine the core and sage datasets
OPT_SETS := $(CURATE)/datasets/filtered-opt.json $(CURATE)/sage/filtered-opt.json
TD_SETS := $(CURATE)/datasets/filtered-td.json $(CURATE)/sage/filtered-td.json $(CURATE)/datasets/filtered-supp-td.json

$(CURATE)/datasets/combined-opt.json: $(OPT_SETS) $(CURATE)/combine.py
	cd $(CURATE) ;					\
	python combine.py combine-opt			\
	 $(addprefix --input-datasets ../,$(OPT_SETS))	\
	--output-dataset ../$@

$(CURATE)/datasets/combined-td.json: $(TD_SETS) $(CURATE)/combine.py
	cd $(CURATE) ;					\
	python combine.py combine-td			\
	$(addprefix --input-datasets ../,$(TD_SETS))	\
	--output-dataset ../$@

## step 2e select parameters from the filtered, combined datasets

### step 2e.i select torsions

DEPS := $(CURATE)/datasets/combined-td.json $(CURATE)/select_parameters.py $(INITIAL_FF)

$(CURATE)/output/td-smirks.json: $(DEPS)
	cd $(CURATE) ;				\
	python select_parameters.py select-td	\
	--dataset datasets/combined-td.json	\
	--forcefield ../$(INITIAL_FF)		\
	--output-smirks output/td-smirks.json	\
	--ring-torsions explicit_ring_torsions.dat

### step 2e.ii select opt
DEPS := $(CURATE)/datasets/combined-opt.json $(CURATE)/select_parameters.py $(INITIAL_FF)

$(CURATE)/output/opt-smirks.json: $(DEPS)
	cd $(CURATE) ;				\
	python select_parameters.py select-opt	\
	--dataset datasets/combined-opt.json	\
	--forcefield ../$(INITIAL_FF)		\
	--output-smirks output/opt-smirks.json



# step 3 - update the initial force field parameters with the modified seminario
# method

$(MSM_FF): $(INITIAL_FF) $(CURATE)/datasets/combined-opt.json $(MSM)/create-msm-ff.py
	cd $(MSM) ; \
	python create-msm-ff.py                                                \
	--initial-force-field       ../$(INITIAL_FF)                           \
	--optimization-dataset      ../$(CURATE)/datasets/combined-opt.json    \
	--working-directory         working-directory                          \
	--output                    ../$@

# step 4 - generate ForceBalance inputs

DEPS := $(FIT)/smiles-to-exclude.dat $(FIT)/smarts-to-exclude.dat		\
	$(CURATE)/datasets/combined-opt.json					\
	$(CURATE)/datasets/combined-td.json $(CURATE)/output/opt-smirks.json	\
	$(CURATE)/output/td-smirks.json $(MSM_FF) $(FIT)/create-fb-inputs.py

$(FIT)/ready: $(DEPS)
	-rm -r $(FIT)/fb-fit/targets
	mkdir -p $(FIT)/fb-fit/targets
	cd $(FIT) ;								\
	python create-fb-inputs.py                                              \
	--tag                       "fb-fit"                                    \
	--optimization-dataset      ../$(CURATE)/datasets/combined-opt.json     \
	--torsion-dataset           ../$(CURATE)/datasets/combined-td.json      \
	--valence-to-optimize       ../$(CURATE)/output/opt-smirks.json		\
	--torsions-to-optimize      ../$(CURATE)/output/td-smirks.json		\
	--forcefield                ../$(MSM_FF)                                \
	--smiles-to-exclude         smiles-to-exclude.dat                       \
	--smarts-to-exclude         smarts-to-exclude.dat                       \
	--max-iterations            100                                         \
	--port                      55387                                       \
	--output-directory          "output"                                    \
	--verbose
	date > $@

step4: $(FIT)/ready

# step 5 - pack up generated files for running on HPC3

# step 5a - pack up the targets directory
$(FIT)/fb-fit/targets.tar.gz: $(FIT)/ready
	-rm $@
	cd $(FIT)/fb-fit ; tar cfz targets.tar.gz targets

# step 5b - zip everything up
TORS_DEPS := $(addprefix $(FIT)/,$(addprefix					\
				        fb-fit/,forcefield/force-field.offxml	\
				        optimize.in targets.tar.gz))

tors.tar.gz: $(FIT)/ready $(TORS_DEPS)
	-rm $@
	tar cfz $@ $(TORS_DEPS) scripts

$(CURATE)/sage/amber-filtered-opt.json: $(CURATE)/amber-charge-filter.py $(CURATE)/sage/opt.json
	python $< --input $(word 2, $^) --output $@

$(CURATE)/sage/amber-filtered-td.json: $(CURATE)/amber-charge-filter.py $(CURATE)/sage/td.json
	python $< --input $(word 2, $^) --output $@
