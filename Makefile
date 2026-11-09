# Reproduce every processed file, table and figure from the raw data.
# Downloads (src/02-04) are not re-run here: raw files are fixed and listed in data/MANIFEST.csv.
#
#   make all      rebuild everything
#   make clean    remove derived outputs

PY      ?= python
HADUK   := data/raw/haduk/tasmax_hadukgrid_uk_country_day_19310101-20251231.nc
NE      := data/raw/naturalearth/ne_10m_admin_0_map_units.zip
DRAWS   ?= 2000

.PHONY: all clean treatment outcomes analysis figures

all: figures

data/processed/tmax_2026_era5_calibrated.csv: src/05_era5_nations.py $(wildcard data/raw/era5/era5_tmax_uk_*.nc)
	$(PY) src/05_era5_nations.py data/raw/era5 $(NE) $(HADUK) config/population_weights.yaml data/processed

treatment: data/processed/tmax_2026_era5_calibrated.csv
	$(PY) src/01_heat_treatment.py $(HADUK) config/population_weights.yaml data/processed data/processed/tmax_2026_era5_calibrated.csv

outcomes:
	$(PY) src/06_build_outcomes.py data data/processed

analysis: treatment outcomes
	$(PY) src/07_analysis.py data/processed config results $(DRAWS)
	$(PY) src/08_size_check.py data/processed results 500

figures: analysis
	$(PY) src/09_figures.py data/processed results figures

clean:
	rm -f data/processed/*.csv results/*.csv results/*.json figures/*.png
