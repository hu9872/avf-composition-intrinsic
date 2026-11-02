# Reproduction driver. Run `make all` after `conda activate avf-composition-intrinsic`
# and `python code/98_relocate_paths.py --to "$PWD"`.

PY := python
CODE := code

.PHONY: all data analysis robustness figures tables verify overlap clean help

help:
	@echo "make data        - download public GEO data (~23 GB) and build AnnData objects"
	@echo "make analysis    - QC, annotation, composition, trajectories, communication, ECM"
	@echo "make robustness  - power/TOST, threshold and gating sensitivity, bulk signature check"
	@echo "make figures     - rebuild main and supplementary figures"
	@echo "make tables      - rebuild supplementary tables and Data S1"
	@echo "make verify      - assert every number in the manuscript against results/"
	@echo "make overlap     - text-overlap pre-screen vs the two source papers"
	@echo "make all         - data + analysis + robustness + figures + tables"

data:
	$(PY) $(CODE)/00_download.py
	$(PY) $(CODE)/00b_download_mt.py
	$(PY) $(CODE)/01_prepare.py

analysis:
	$(PY) $(CODE)/02_qc.py
	$(PY) $(CODE)/03_annotate.py
	$(PY) $(CODE)/04_fibro.py
	$(PY) $(CODE)/05_trajectory.py
	$(PY) $(CODE)/05b_composition_pseudobulk.py
	$(PY) $(CODE)/05c_compartment_pseudobulk.py
	$(PY) $(CODE)/05d_myeloid.py
	$(PY) $(CODE)/05e_cdc2.py
	$(PY) $(CODE)/05f_decomposition.py
	$(PY) $(CODE)/05g_gene_origin_all.py
	$(PY) $(CODE)/06_communication.py
	$(PY) $(CODE)/06c_pathways.py
	$(PY) $(CODE)/07_ecm.py
	$(PY) $(CODE)/08_bulk.py
	$(PY) $(CODE)/08b_vein_to_avf.py

robustness:
	$(PY) $(CODE)/30_power_tost.py
	$(PY) $(CODE)/31_rule_floor_sensitivity.py
	$(PY) $(CODE)/32_gating_perturbation.py
	$(PY) $(CODE)/33_deconvolution.py
	$(PY) $(CODE)/34_compartment_signatures.py
	$(PY) $(CODE)/35_extra_checks.py

figures:
	$(PY) $(CODE)/09_figures.py
	$(PY) $(CODE)/09d_figures_myeloid.py
	$(PY) $(CODE)/18_figures_final.py
	$(PY) $(CODE)/40_graphical_abstract.py

tables:
	$(PY) $(CODE)/20_supp_tables.py
	$(PY) $(CODE)/44_make_supp_tables.py

verify:
	$(PY) $(CODE)/90_verify_numbers.py

overlap:
	$(PY) $(CODE)/95_text_overlap_check.py

all: data analysis robustness figures tables verify

clean:
	find . -name "__pycache__" -type d -prune -exec rm -rf {} +
