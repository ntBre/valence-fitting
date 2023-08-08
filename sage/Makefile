# Refit the original Sage force field on the new data set

GENERATE := 01_generate-forcefield
CURATE := 02_curate-data
MSM := 03_generate-initial-ff
FIT := 04_fit-forcefield


# step 1 - generate force field - is nearly done, simply load the force field
# and dump it back out
initial.offxml: initial.py
	python initial.py

INITIAL_FF := initial.offxml


# step 2 - curate data set - is the same process, but with a different
# --initial-forcefield argument to select_parameters, so call the original
# makefile up to that point

../$(CURATE)/datasets/filtered-td.json:
	cd .. ; \
	make $(CURATE)/datasets/filtered-td.json

../$(CURATE)/datasets/filtered-opt.json:
	cd .. ; \
	make $(CURATE)/datasets/filtered-opt.json

DEPS := ../$(CURATE)/datasets/filtered-td.json		\
	../$(CURATE)/select_parameters.py $(INITIAL_FF)

td-smirks.json: $(DEPS)
	python ../$(CURATE)/select_parameters.py select-td	\
	--dataset ../$(CURATE)/datasets/filtered-td.json	\
	--forcefield $(INITIAL_FF)				\
	--output-smirks $@					\
	--ring-torsions ../$(CURATE)/explicit_ring_torsions.dat

### step 2d.ii select opt
DEPS := ../$(CURATE)/datasets/filtered-opt.json		\
	../$(CURATE)/select_parameters.py $(INITIAL_FF)

opt-smirks.json: $(DEPS)
	python ../$(CURATE)/select_parameters.py select-opt	\
	--dataset ../$(CURATE)/datasets/filtered-opt.json	\
	--forcefield $(INITIAL_FF)				\
	--output-smirks $@


# step 3 - generate initial parameters using MSM guess

MSM_FF := msm.offxml

$(MSM_FF): $(INITIAL_FF) ../$(CURATE)/datasets/filtered-opt.json ../$(MSM)/create-msm-ff.py
	python ../$(MSM)/create-msm-ff.py                                       \
	--initial-force-field       $(INITIAL_FF)				\
	--optimization-dataset      ../$(CURATE)/output/filtered-opt.json	\
	--working-directory         working-directory				\
	--output                    $@


# step 4 - generate ForceBalance inputs

DEPS := ../$(FIT)/smiles-to-exclude.dat ../$(FIT)/smarts-to-exclude.dat		\
	../$(CURATE)/datasets/filtered-opt.json					\
	../$(CURATE)/datasets/filtered-td.json opt-smirks.json td-smirks.json	\
	$(MSM_FF) ../$(FIT)/create-fb-inputs.py

step4: $(DEPS)
	rm -r fb-fit/targets
	mkdir -p fb-fit/targets
	python ../$(FIT)/create-fb-inputs.py                                    \
	--tag                       "fb-fit"                                    \
	--optimization-dataset      ../$(CURATE)/datasets/filtered-opt.json     \
	--torsion-dataset           ../$(CURATE)/datasets/filtered-td.json      \
	--valence-to-optimize       opt-smirks.json				\
	--torsions-to-optimize      td-smirks.json				\
	--forcefield                $(MSM_FF)					\
	--smiles-to-exclude         ../$(FIT)/smiles-to-exclude.dat             \
	--smarts-to-exclude         ../$(FIT)/smarts-to-exclude.dat             \
	--max-iterations            1                                           \
	--port                      55387                                       \
	--output-directory          "output"                                    \
	--verbose
	date > $@

fb-fit/targets.tar.gz: $(wildcard fb-fit/targets/*/*)
	rm $@
	cd fb-fit ; tar cfz targets.tar.gz targets

DEPS := $(addprefix fb-fit/,forcefield/force-field.offxml optimize.in	\
				        targets.tar.gz)

sage.tar.gz: step4 $(DEPS)
	rm $@
	cp -r ../scripts .
	tar cfz $@ $(DEPS) scripts
	rm -r scripts
