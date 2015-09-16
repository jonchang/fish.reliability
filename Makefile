PCS = PC1.txt PC2.txt PC3.txt PC4.txt PC5.txt PC6.txt
DATA = data/fish_families.rda data/fish_reliability.rda

all: data bamm

.PHONY: clean

clean:
	rm $(DATA)

data: $(DATA)

data/%.rda: data-raw/%.R
	cd '$(<D)' && Rscript '$(<F)'

bamm: data-raw/PC1.txt

data-raw/PC1.txt: data-raw/bamm_pca.R data/fish_families.rda
	cd data-raw && Rscript bamm_pca.R
