# Refit the original Sage force field on the new data set

GENERATE := 01_generate-forcefield
CURATE := 02_curate-data
MSM := 03_generate-initial-ff
FIT := 04_fit-forcefield

.PHONY: step3

# step 1 - generate force field - is nearly done, simply load the force field
# and dump it back out
initial.offxml: initial.py
	python initial.py

INITIAL_FF := initial.offxml

opt := ../$(CURATE)/datasets/combined-opt.json
td := ../$(CURATE)/datasets/combined-td.json

# step 2 - curate data set - is the same process, but with a different
# --initial-forcefield argument to select_parameters, so call the original
# makefile up to that point

$(opt):
	make -C .. $(CURATE)/datasets/combined-td.json

$(td):
	make -C .. $(CURATE)/datasets/combined-opt.json

td-smirks.json: $(td) ../$(CURATE)/select_parameters.py $(INITIAL_FF)
	python ../$(CURATE)/select_parameters.py select-td	\
	--dataset $(td)						\
	--forcefield $(INITIAL_FF)				\
	--output-smirks $@					\
	--ring-torsions ../$(CURATE)/explicit_ring_torsions.dat

### step 2d.ii select opt
opt-smirks.json: $(opt) ../$(CURATE)/select_parameters.py $(INITIAL_FF)
	python ../$(CURATE)/select_parameters.py select-opt	\
	--dataset $(opt)					\
	--forcefield $(INITIAL_FF)				\
	--output-smirks $@

# step 3 - generate initial parameters using MSM guess

MSM_FF := msm.offxml

$(MSM_FF): $(INITIAL_FF) $(opt) ../$(MSM)/create-msm-ff.py
	python ../$(MSM)/create-msm-ff.py               \
	--initial-force-field       $(INITIAL_FF)	\
	--optimization-dataset      $(opt)		\
	--working-directory         working-directory	\
	--output                    $@

step3: $(MSM_FF)

# step 4 - generate ForceBalance inputs
DEPS := ../$(FIT)/smiles-to-exclude.dat ../$(FIT)/smarts-to-exclude.dat $(opt)	\
	$(td) opt-smirks.json td-smirks.json $(MSM_FF)				\
	../$(FIT)/create-fb-inputs.py

fb_inputs := $(addprefix fb-fit/,forcefield/force-field.offxml optimize.in targets)

$(fb_inputs): $(DEPS)
	rm -rf fb-fit/targets
	mkdir -p fb-fit/targets
	python ../$(FIT)/create-fb-inputs.py                                    \
	--tag                       "fb-fit"                                    \
	--optimization-dataset      $(opt)     \
	--torsion-dataset           $(td)      \
	--valence-to-optimize       opt-smirks.json				\
	--torsions-to-optimize      td-smirks.json				\
	--forcefield                $(MSM_FF)					\
	--smiles-to-exclude         ../$(FIT)/smiles-to-exclude.dat             \
	--smarts-to-exclude         ../$(FIT)/smarts-to-exclude.dat             \
	--max-iterations            100                                         \
	--port                      55387                                       \
	--output-directory          "output"                                    \
	--verbose

fb-fit/targets.tar.gz: $(wildcard fb-fit/targets/*/*)
	rm -f $@
	cd fb-fit ; tar cfz targets.tar.gz targets

fb_inputs := $(addprefix fb-fit/,forcefield/force-field.offxml optimize.in targets.tar.gz)
sage.tar.gz: $(fb_inputs)
	rm -f $@
	cp -r ../scripts .
	tar cfz $@ $^ scripts
	rm -r scripts
