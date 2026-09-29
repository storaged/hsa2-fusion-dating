# Quick targets. All paths are relative to the repository root.
PYTHON ?= python
FIG = workflow/09_figures
OUT = results/figures

.PHONY: help figures test install clean-figures

help:
	@echo "make figures   regenerate all figures, supplementary tables and numbers.tex from the shipped results (<1 min)"
	@echo "make test      smoke tests: model identities and small exact/simulation runs (~1 min)"
	@echo "make install   copy the workflow into \$$HSA2_ROOT/scripts for a full re-run (see REPRODUCE.md)"

figures:
	mkdir -p $(OUT)
	$(PYTHON) -W ignore $(FIG)/F_schematic.py $(OUT)
	$(PYTHON) -W ignore $(FIG)/F_make_figures.py
	$(PYTHON) -W ignore $(FIG)/F_supp_figures.py
	$(PYTHON) -W ignore $(FIG)/F_supp_tables.py
	$(PYTHON) $(FIG)/gen_numbers.py $(OUT)/stats.json $(OUT)/numbers.tex
	@echo "figures written to $(OUT)/"

test:
	bash tests/smoke.sh

install:
	utils/install_scripts.sh

clean-figures:
	rm -f $(OUT)/Fig*.pdf $(OUT)/Fig*.png $(OUT)/supp_tables.tex $(OUT)/numbers.tex $(OUT)/Supplementary_Data_1.xlsx
