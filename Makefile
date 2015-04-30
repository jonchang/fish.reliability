data: data/fish_families.rda data/fish_reliability.rda

data/%.rda: data-raw/%.R
	cd '$(<D)' && Rscript '$(<F)'
