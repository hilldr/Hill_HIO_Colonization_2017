## LaTeX Makefile
## define shorthand file names for text
TEXT=./src/elife_format
TEXTNF=./src/eLife_nofigures
TEXTBR=./src/complete_paper
FINALTEXT=Hill_eLife
BIORXIV=Hill_bioRxiv

## output to eLife PDF
pdf: $(FINALTEXT).pdf
$(FINALTEXT).pdf: $(TEXT).tex ./src/bibliography.bib \
	./figures/figure1/figure1_multipanel.pdf \
	./figures/figure2/figure2_multipanel.pdf \
	./figures/figure3/figure3_multipanel.pdf \
	./figures/figure4/figure4_multipanel.pdf \
	./figures/figure5/figure5_multipanel.pdf \
	./figures/figure6/figure6_multipanel.pdf \
	./figures/figure7/figure7_multipanel.pdf \
	./figures/figure8/figure8_multipanel.pdf		
	pdflatex -output-directory src $(TEXT)
	pdflatex -output-directory src $(TEXT)
	cp ./src/bibliography.bib ./
	cp ./src/elife.cls ./
	cp ./src/vancouver-elife.bst ./
	bibtex $(TEXT)
	pdflatex -output-directory src $(TEXT)
	pdflatex -output-directory src $(TEXT)
	mv $(TEXT).pdf $(FINALTEXT).pdf
	rm *.bib *.cls *.bst

## output "track changes" pdf
changes: $(FINALTEXT)_changes.pdf
$(FINALTEXT)_changes.pdf: $(TEXT).tex ./src/bibliography.bib \
	./figures/figure1/figure1_multipanel.pdf \
	./figures/figure2/figure2_multipanel.pdf \
	./figures/figure3/figure3_multipanel.pdf \
	./figures/figure4/figure4_multipanel.pdf \
	./figures/figure5/figure5_multipanel.pdf \
	./figures/figure6/figure6_multipanel.pdf \
	./figures/figure7/figure7_multipanel.pdf \
	./figures/figure8/figure8_multipanel.pdf \
	./revisions/elife-2nd-submission.tex
	latexdiff --type=UNDERLINE -s COLOR ./revisions/elife-2nd-submission.tex $(TEXT).tex > $(TEXT)_changes.tex
	pdflatex -output-directory src $(TEXT)_changes
	pdflatex -output-directory src $(TEXT)_changes
	cp ./src/bibliography.bib ./
	cp ./src/elife.cls ./
	cp ./src/vancouver-elife.bst ./
	bibtex $(TEXT)_changes
	pdflatex -output-directory src $(TEXT)_changes
	pdflatex -output-directory src $(TEXT)_changes
	mv $(TEXT)_changes.pdf $(FINALTEXT)_changes.pdf
	rm *.bib *.cls *.bst 

## output to bioRxiv PDF
## work in progress, partially working
br-pdf: $(BIORXIV).pdf
$(BIORXIV).pdf: $(TEXTBR).tex ./src/bibliography.bib \
	./figures/figure1/figure1_multipanel.pdf \
	./figures/figure2/figure2_multipanel.pdf \
	./figures/figure3/figure3_multipanel.pdf \
	./figures/figure4/figure4_multipanel.pdf \
	./figures/figure5/figure5_multipanel.pdf \
	./figures/figure6/figure6_multipanel.pdf \
	./figures/figure7/figure7_multipanel.pdf \
	./figures/figure8/figure8_multipanel.pdf		
	pdflatex -output-directory src $(TEXTBR)
	pdflatex -output-directory src $(TEXTBR)
	cp ./src/bibliography.bib ./
	cp ./src/elsarticle.cls ./
	cp ./src/model2-names.bst ./
	bibtex $(TEXTBR)
	pdflatex -output-directory src $(TEXTBR)
	pdflatex -output-directory src $(TEXTBR)
	mv $(TEXTBR).pdf $(BIORXIV).pdf
	rm *.bib *.cls *.bst

## output compressed PDF
DENSITY=50
QUALITY=50
compressed: $(FINALTEXT)_compressed.pdf
$(FINALTEXT)_compressed.pdf: $(TEXT).tex ./src/bibliography.bib \
	./figures/figure1/figure1_multipanel.pdf \
	./figures/figure2/figure2_multipanel.pdf \
	./figures/figure3/figure3_multipanel.pdf \
	./figures/figure4/figure4_multipanel.pdf \
	./figures/figure5/figure5_multipanel.pdf \
	./figures/figure6/figure6_multipanel.pdf \
	./figures/figure7/figure7_multipanel.pdf \
	./figures/figure8/figure8_multipanel.pdf
	convert -density $(DENSITY) ./figures/figure1/figure1_multipanel.pdf -quality $(QUALITY) ./figures/figure1/figure1_multipanel.png
	convert -density $(DENSITY) ./figures/figure2/figure2_multipanel.pdf -quality $(QUALITY) ./figures/figure2/figure2_multipanel.png
	convert -density $(DENSITY) ./figures/figure3/figure3_multipanel.pdf -quality $(QUALITY) ./figures/figure3/figure3_multipanel.png
	convert -density $(DENSITY) ./figures/figure4/figure4_multipanel.pdf -quality $(QUALITY) ./figures/figure4/figure4_multipanel.png
	convert -density $(DENSITY) ./figures/figure5/figure5_multipanel.pdf -quality $(QUALITY) ./figures/figure5/figure5_multipanel.png
	convert -density $(DENSITY) ./figures/figure6/figure6_multipanel.pdf -quality $(QUALITY) ./figures/figure6/figure6_multipanel.png
	convert -density $(DENSITY) ./figures/figure7/figure7_multipanel.pdf -quality $(QUALITY) ./figures/figure7/figure7_multipanel.png
	convert -density $(DENSITY) ./figures/figure8/figure8_multipanel.pdf -quality $(QUALITY) ./figures/figure8/figure8_multipanel.png
	# convert -density $(DENSITY) ./figures/supplemental/sfigure1_multipanel.pdf -quality $(QUALITY) ./figures/supplemental/sfigure1_multipanel.png
	# convert -density $(DENSITY) ./figures/figure1/sfigure1-2_tree.pdf -quality $(QUALITY) ./figures/figure1/sfigure1-2_tree.png
	# convert -density $(DENSITY) ./figures/figure1/sfigure1-3_multipanel.pdf -quality $(QUALITY) ./figures/figure1/sfigure1-3_multipanel.png
	# convert -density $(DENSITY) ./figures/supplemental/sfigure1_supp2.pdf -quality $(QUALITY) ./figures/supplemental/sfigure1_supp2.png
	# convert -density $(DENSITY) ./results/ECOR2HIO_24-96-RNAseq/KEGG-pathview/Fig3s1_hsa04110_cellcycle_48.pdf -quality $(QUALITY) ./results/ECOR2HIO_24-96-RNAseq/KEGG-pathview/Fig3s1_hsa04110_cellcycle_48.png
	# convert -density $(DENSITY) ./figures/supplemental/figure5s1_multipanel.pdf -quality $(QUALITY) ./figures/supplemental/figure5s1_multipanel.png
	# convert -density $(DENSITY) ./figures/supplemental/sfigure4_supp3.pdf -quality $(QUALITY) ./figures/supplemental/sfigure4_supp3.png
	# convert -density $(DENSITY) ./figures/supplemental/DEFB4B-expression.pdf -quality $(QUALITY) ./figures/supplemental/DEFB4B-expression.png
	# convert -density $(DENSITY) ./figures/figure6/figure6_supplement2.pdf -quality $(QUALITY) ./figures/figure6/figure6_supplement2.png
	# convert -density $(DENSITY) ./figures/figure7/Supplemental_Figure4_Muc2-NFkB.pdf -quality $(QUALITY) ./figures/figure7/Supplemental_Figure4_Muc2-NFkB.png
	# convert -density $(DENSITY) ./figures/supplemental/sfigure4_multipanel.pdf -quality $(QUALITY) ./figures/supplemental/sfigure4_multipanel.png
	cp $(TEXT).tex $(TEXT)_compressed.tex
	sed -i 's/\.pdf/\.png/g' $(TEXT)_compressed.tex # use png versions of figures
	pdflatex -output-directory src $(TEXT)_compressed
	pdflatex -output-directory src $(TEXT)_compressed
	cp ./src/bibliography.bib ./
	cp ./src/elife.cls ./
	cp ./src/vancouver-elife.bst ./
	bibtex $(TEXT)_compressed
	pdflatex -output-directory src $(TEXT)_compressed
	pdflatex -output-directory src $(TEXT)_compressed
	mv $(TEXT)_compressed.pdf $(FINALTEXT)_compressed.pdf
	rm *.bib *.cls *.bst

## output compressed PDF
DENSITY=50
QUALITY=50
compressed-full: $(FINALTEXT)_compressed_full.pdf
$(FINALTEXT)_compressed_full.pdf: $(TEXT).tex ./src/bibliography.bib \
	./figures/figure1/figure1_multipanel.png \
	./figures/figure2/figure2_multipanel.png \
	./figures/figure3/figure3_multipanel.png \
	./figures/figure4/figure4_multipanel.png \
	./figures/figure5/figure5_multipanel.png \
	./figures/figure6/figure6_multipanel.png \
	./figures/figure7/figure7_multipanel.png \
	./figures/figure8/figure8_multipanel.png \
	./figures/supplemental/sfigure1_multipanel.png \
	./figures/figure1/sfigure1-2_tree.png \
	./figures/figure1/sfigure1-3_multipanel.png \
	./results/ECOR2HIO_24-96-RNAseq/KEGG-pathview/Fig3s1_hsa04110_cellcycle_48.png \
	./figures/supplemental/figure5s1_multipanel.png \
	./figures/supplemental/sfigure4_supp3.png \
	./figures/supplemental/DEFB4B-expression.png \
	./figures/figure6/figure6_supplement2.png \
	./figures/figure7/Supplemental_Figure4_Muc2-NFkB.png \
	./figures/supplemental/sfigure4_multipanel.png 
	cp $(TEXT).tex $(TEXT)_compressed_full.tex
	sed -i 's/\.pdf/.png/g' $(TEXT)_compressed_full.tex # use png versions of figures
	sed -i "`wc -l < $(TEXT)_compressed_full.tex`i\\\includepdf{./supplements/figure1_supplement1.pdf}\\" $(TEXT)_compressed_full.tex 
	sed -i "`wc -l < $(TEXT)_compressed_full.tex`i\\\includepdf{./supplements/figure1_supplement2.pdf}\\" $(TEXT)_compressed_full.tex
	sed -i "`wc -l < $(TEXT)_compressed_full.tex`i\\\includepdf{./supplements/figure1_supplement3.pdf}\\" $(TEXT)_compressed_full.tex  
	sed -i "`wc -l < $(TEXT)_compressed_full.tex`i\\\includepdf{./supplements/figure3_supplement1.pdf}\\" $(TEXT)_compressed_full.tex  
	sed -i "`wc -l < $(TEXT)_compressed_full.tex`i\\\includepdf{./supplements/figure5_supplement1.pdf}\\" $(TEXT)_compressed_full.tex
	sed -i "`wc -l < $(TEXT)_compressed_full.tex`i\\\includepdf{./supplements/figure5_supplement2.pdf}\\" $(TEXT)_compressed_full.tex
	sed -i "`wc -l < $(TEXT)_compressed_full.tex`i\\\includepdf{./supplements/figure6_supplement1.pdf}\\" $(TEXT)_compressed_full.tex
	sed -i "`wc -l < $(TEXT)_compressed_full.tex`i\\\includepdf{./supplements/figure6_supplement2.pdf}\\" $(TEXT)_compressed_full.tex
	sed -i "`wc -l < $(TEXT)_compressed_full.tex`i\\\includepdf{./supplements/figure7_supplement1.pdf}\\" $(TEXT)_compressed_full.tex
	sed -i "`wc -l < $(TEXT)_compressed_full.tex`i\\\includepdf{./supplements/figure8_supplement1.pdf}\\" $(TEXT)_compressed_full.tex
	pdflatex -output-directory src $(TEXT)_compressed_full
	pdflatex -output-directory src $(TEXT)_compressed_full
	cp ./src/bibliography.bib ./
	cp ./src/elife.cls ./
	cp ./src/vancouver-elife.bst ./
	bibtex $(TEXT)_compressed_full
	pdflatex -output-directory src $(TEXT)_compressed_full
	pdflatex -output-directory src $(TEXT)_compressed_full
	mv $(TEXT)_compressed_full.pdf $(FINALTEXT)_compressed_full.pdf
	rm *.bib *.cls *.bst

## output to eLife PDF
pdf-nofigures: $(FINALTEXT)_noFigures.pdf
$(FINALTEXT)_noFigures.pdf: $(TEXTNF).tex ./src/bibliography.bib 
	pdflatex -output-directory src $(TEXTNF)
	pdflatex -output-directory src $(TEXTNF)
	cp ./src/bibliography.bib ./
	cp ./src/elife.cls ./
	cp ./src/vancouver-elife.bst ./
	bibtex $(TEXTNF)
	pdflatex -output-directory src $(TEXTNF)
	pdflatex -output-directory src $(TEXTNF)
	mv $(TEXTNF).pdf $(FINALTEXT)_noFigures.pdf
	rm *.bib *.cls *.bst

## make PNG versions of figures
./figures/figure1/figure1_multipanel.png: ./figures/figure1/figure1_multipanel.pdf
	convert -density $(DENSITY) ./figures/figure1/figure1_multipanel.pdf -quality $(QUALITY) ./figures/figure1/figure1_multipanel.png
./figures/figure2/figure2_multipanel.png: ./figures/figure2/figure2_multipanel.pdf
	convert -density $(DENSITY) ./figures/figure2/figure2_multipanel.pdf -quality $(QUALITY) ./figures/figure2/figure2_multipanel.png
./figures/figure3/figure3_multipanel.png: ./figures/figure3/figure3_multipanel.pdf
	convert -density $(DENSITY) ./figures/figure3/figure3_multipanel.pdf -quality $(QUALITY) ./figures/figure3/figure3_multipanel.png
./figures/figure4/figure4_multipanel.png: ./figures/figure4/figure4_multipanel.pdf
	convert -density $(DENSITY) ./figures/figure4/figure4_multipanel.pdf -quality $(QUALITY) ./figures/figure4/figure4_multipanel.png
./figures/figure5/figure5_multipanel.png: ./figures/figure5/figure5_multipanel.pdf
	convert -density $(DENSITY) ./figures/figure5/figure5_multipanel.pdf -quality $(QUALITY) ./figures/figure5/figure5_multipanel.png
./figures/figure6/figure6_multipanel.png: ./figures/figure6/figure6_multipanel.pdf
	convert -density $(DENSITY) ./figures/figure6/figure6_multipanel.pdf -quality $(QUALITY) ./figures/figure6/figure6_multipanel.png
./figures/figure7/figure7_multipanel.png: ./figures/figure7/figure7_multipanel.pdf
	convert -density $(DENSITY) ./figures/figure7/figure7_multipanel.pdf -quality $(QUALITY) ./figures/figure7/figure7_multipanel.png
./figures/figure8/figure8_multipanel.png: ./figures/figure8/figure8_multipanel.pdf
	convert -density $(DENSITY) ./figures/figure8/figure8_multipanel.pdf -quality $(QUALITY) ./figures/figure8/figure8_multipanel.png
./figures/supplemental/sfigure1_multipanel.png: ./figures/supplemental/sfigure1_multipanel.pdf
	convert -density $(DENSITY) ./figures/supplemental/sfigure1_multipanel.pdf -quality $(QUALITY) ./figures/supplemental/sfigure1_multipanel.png
./figures/figure1/sfigure1-2_tree.png: ./figures/figure1/sfigure1-2_tree.pdf
	convert -density $(DENSITY) ./figures/figure1/sfigure1-2_tree.pdf -quality $(QUALITY) ./figures/figure1/sfigure1-2_tree.png
./figures/figure1/sfigure1-3_multipanel.png: ./figures/figure1/sfigure1-3_multipanel.pdf
	convert -density $(DENSITY) ./figures/figure1/sfigure1-3_multipanel.pdf -quality $(QUALITY) ./figures/figure1/sfigure1-3_multipanel.png
./figures/supplemental/sfigure1_supp2.png: ./figures/supplemental/sfigure1_supp2.pdf
	convert -density $(DENSITY) ./figures/supplemental/sfigure1_supp2.pdf -quality $(QUALITY) ./figures/supplemental/sfigure1_supp2.png
./results/ECOR2HIO_24-96-RNAseq/KEGG-pathview/Fig3s1_hsa04110_cellcycle_48.png: ./results/ECOR2HIO_24-96-RNAseq/KEGG-pathview/Fig3s1_hsa04110_cellcycle_48.pdf
	convert -density $(DENSITY) ./results/ECOR2HIO_24-96-RNAseq/KEGG-pathview/Fig3s1_hsa04110_cellcycle_48.pdf -quality $(QUALITY) ./results/ECOR2HIO_24-96-RNAseq/KEGG-pathview/Fig3s1_hsa04110_cellcycle_48.png
./figures/supplemental/figure5s1_multipanel.png: ./figures/supplemental/figure5s1_multipanel.pdf
	convert -density $(DENSITY) ./figures/supplemental/figure5s1_multipanel.pdf -quality $(QUALITY) ./figures/supplemental/figure5s1_multipanel.png
./figures/supplemental/sfigure4_supp3.png: ./figures/supplemental/sfigure4_supp3.pdf
	convert -density $(DENSITY) ./figures/supplemental/sfigure4_supp3.pdf -quality $(QUALITY) ./figures/supplemental/sfigure4_supp3.png
./figures/supplemental/DEFB4B-expression.png: ./figures/supplemental/DEFB4B-expression.pdf
	convert -density $(DENSITY) ./figures/supplemental/DEFB4B-expression.pdf -quality $(QUALITY) ./figures/supplemental/DEFB4B-expression.png
./figures/figure6/figure6_supplement2.png: ./figures/figure6/figure6_supplement2.pdf
	convert -density $(DENSITY) ./figures/figure6/figure6_supplement2.pdf -quality $(QUALITY) ./figures/figure6/figure6_supplement2.png
./figures/figure7/Supplemental_Figure4_Muc2-NFkB.png: ./figures/figure7/Supplemental_Figure4_Muc2-NFkB.pdf
	convert -density $(DENSITY) ./figures/figure7/Supplemental_Figure4_Muc2-NFkB.pdf -quality $(QUALITY) ./figures/figure7/Supplemental_Figure4_Muc2-NFkB.png
./figures/supplemental/sfigure4_multipanel.png: ./figures/supplemental/sfigure4_multipanel.pdf
	convert -density $(DENSITY) ./figures/supplemental/sfigure4_multipanel.pdf -quality $(QUALITY) ./figures/supplemental/sfigure4_multipanel.png


## supplemental figures
supplements: ./supplements/figure1_supplement1.pdf \
	./supplements/figure1_supplement2.pdf \
	./supplements/figure1_supplement3.pdf \
	./supplements/figure3_supplement1.pdf \
	./supplements/figure5_supplement1.pdf \
	./supplements/figure5_supplement2.pdf \
	./supplements/figure6_supplement1.pdf \
	./supplements/figure6_supplement2.pdf \
	./supplements/figure7_supplement1.pdf \
	./supplements/figure8_supplement1.pdf 

./supplements/figure1_supplement1.pdf: ./src/figure1_supplement1.tex
	pdflatex -output-directory src ./src/figure1_supplement1
	pdflatex -output-directory src ./src/figure1_supplement1
	cp ./src/bibliography.bib ./
	cp ./src/elsarticle.cls ./
	cp ./src/vancouver-elife.bst ./
	mv ./src/figure1_supplement1.pdf ./supplements/figure1_supplement1.pdf
	rm *.bib *.cls *.bst

./supplements/figure1_supplement2.pdf: ./src/figure1_supplement2.tex
	pdflatex -output-directory src ./src/figure1_supplement2
	pdflatex -output-directory src ./src/figure1_supplement2
	cp ./src/bibliography.bib ./
	cp ./src/elsarticle.cls ./
	cp ./src/vancouver-elife.bst ./
	bibtex ./src/figure1_supplement2
	pdflatex -output-directory src ./src/figure1_supplement2
	pdflatex -output-directory src ./src/figure1_supplement2
	mv ./src/figure1_supplement2.pdf ./supplements/figure1_supplement2.pdf
	rm *.bib *.cls *.bst

./supplements/figure1_supplement3.pdf: ./src/figure1_supplement3.tex
	pdflatex -output-directory src ./src/figure1_supplement3
	pdflatex -output-directory src ./src/figure1_supplement3
	cp ./src/bibliography.bib ./
	cp ./src/elsarticle.cls ./
	cp ./src/vancouver-elife.bst ./
	mv ./src/figure1_supplement3.pdf ./supplements/figure1_supplement3.pdf
	rm *.bib *.cls *.bst

./supplements/figure3_supplement1.pdf: ./src/figure3_supplement1.tex
	pdflatex -output-directory src ./src/figure3_supplement1
	pdflatex -output-directory src ./src/figure3_supplement1
	cp ./src/bibliography.bib ./
	cp ./src/elsarticle.cls ./
	cp ./src/vancouver-elife.bst ./
	mv ./src/figure3_supplement1.pdf ./supplements/figure3_supplement1.pdf
	rm *.bib *.cls *.bst

./supplements/figure5_supplement1.pdf: ./src/figure5_supplement1.tex
	pdflatex -output-directory src ./src/figure5_supplement1
	pdflatex -output-directory src ./src/figure5_supplement1
	cp ./src/bibliography.bib ./
	cp ./src/elsarticle.cls ./
	cp ./src/vancouver-elife.bst ./
	mv ./src/figure5_supplement1.pdf ./supplements/figure5_supplement1.pdf
	rm *.bib *.cls *.bst

./supplements/figure5_supplement2.pdf: ./src/figure5_supplement2.tex
	pdflatex -output-directory src ./src/figure5_supplement2
	pdflatex -output-directory src ./src/figure5_supplement2
	cp ./src/bibliography.bib ./
	cp ./src/elsarticle.cls ./
	cp ./src/vancouver-elife.bst ./
	mv ./src/figure5_supplement2.pdf ./supplements/figure5_supplement2.pdf
	rm *.bib *.cls *.bst

./supplements/figure6_supplement1.pdf: ./src/figure6_supplement1.tex
	pdflatex -output-directory src ./src/figure6_supplement1
	pdflatex -output-directory src ./src/figure6_supplement1
	cp ./src/bibliography.bib ./
	cp ./src/elsarticle.cls ./
	cp ./src/vancouver-elife.bst ./
	mv ./src/figure6_supplement1.pdf ./supplements/figure6_supplement1.pdf
	rm *.bib *.cls *.bst

./supplements/figure6_supplement2.pdf: ./src/figure6_supplement2.tex
	pdflatex -output-directory src ./src/figure6_supplement2
	pdflatex -output-directory src ./src/figure6_supplement2
	cp ./src/bibliography.bib ./
	cp ./src/elsarticle.cls ./
	cp ./src/vancouver-elife.bst ./
	mv ./src/figure6_supplement2.pdf ./supplements/figure6_supplement2.pdf
	rm *.bib *.cls *.bst

./supplements/figure7_supplement1.pdf: ./src/figure7_supplement1.tex
	pdflatex -output-directory src ./src/figure7_supplement1
	pdflatex -output-directory src ./src/figure7_supplement1
	cp ./src/bibliography.bib ./
	cp ./src/elsarticle.cls ./
	cp ./src/vancouver-elife.bst ./
	mv ./src/figure7_supplement1.pdf ./supplements/figure7_supplement1.pdf
	rm *.bib *.cls *.bst

./supplements/figure8_supplement1.pdf: ./src/figure8_supplement1.tex
	pdflatex -output-directory src ./src/figure8_supplement1
	pdflatex -output-directory src ./src/figure8_supplement1
	cp ./src/bibliography.bib ./
	cp ./src/elsarticle.cls ./
	cp ./src/vancouver-elife.bst ./
	mv ./src/figure8_supplement1.pdf ./supplements/figure8_supplement1.pdf
	rm *.bib *.cls *.bst

## output oly
figures: $(FINALTEXT)_figures.pdf
$(FINALTEXT)_figures.pdf : ./src/figures_only.tex
	pdflatex -output-directory src ./src/figures_only
	mv ./src/figures_only.pdf $(FINALTEXT)figures.pdf

## output to DOCX
docx: $(FINALTEXT).docx
$(FINALTEXT).docx: $(TEXT).tex
	cp $(TEXT).tex $(TEXT)_docx_reformat.tex
	sed -i 's/pdf/png/g' $(TEXT)_docx_reformat.tex # use png versions of figures
	sed -i 's/\\(\\kappa\\)/κ/g' $(TEXT)_docx_reformat.tex
	sed -i 's/\\(\\beta\\)/β/g' $(TEXT)_docx_reformat.tex
	sed -i 's/\\(\\alpha\\)/α/g' $(TEXT)_docx_reformat.tex
	sed -i 's/\\(\\mu\\)/μ/g' $(TEXT)_docx_reformat.tex
	sed -i 's/\\(\\gamma\\)/γ/g' $(TEXT)_docx_reformat.tex
	sed -i 's/{\"i}/ï/g' $(TEXT)_docx_reformat.tex
	sed -i 's/\\pm/±/g' $(TEXT)_docx_reformat.tex
	sed -i 's/\num{//g' $(TEXT)_docx_reformat.tex
	sed -i 's/\\(_{\\text{2}}\\)/₂/g' $(TEXT)_docx_reformat.tex
	sed -i 's/\\(^{\\text{2}}\\)/²/g' $(TEXT)_docx_reformat.tex
	pandoc --bibliography=./src/bibliography.bib --filter pandoc-citeproc  $(TEXT)_docx_reformat.tex -o $(FINALTEXT).docx

.PHONY: clean
clean:
	cd src && rm *.aux *.blg *.out *.bbl *.log

## Generate Figure 1 multipanel figure
figure1: ./figures/figure1/figure1_multipanel.pdf
./figures/figure1/figure1_multipanel.pdf : ./src/figure_Rscripts/figure1.R \
	./data/figure1/010716_01_R3D.csv \
	./data/figure1/010716/010716_01_R3D_w594_t01.png \
	./figures/figure1/figure1b.png \
	./data/figure1/ECOR2growth_fig1.csv \
	./data/figure1/sample_table_fig1.csv \
	./data/figure1/ECOR2growth_fig1_timecourse.csv \
	./data/figure1/161206_survival/survival_and_ELISA.csv
	R -e "setwd('./src/'); source('figure_Rscripts/figure1.R')"

## Generate gene counts and DE from kallisto output for Figure 2
./results/ECOR2HIO_24-96-RNAseq/complete-dataset_DESeq2-normalized-counts.csv ./results/ECOR2HIO_24-96-RNAseq/ECOR2_over_PBS_24hr.csv ./results/ECOR2HIO_24-96-RNAseq/ECOR2_over_PBS_48hr.csv ./results/ECOR2HIO_24-96-RNAseq/ECOR2_over_PBS_96hr.csv : ./src/figure_Rscripts/figure3-kallisto-post.R \
	data/RNA-seq/kallisto-2014/Sample_34966_R1.fastq/abundance.tsv \
	data/RNA-seq/kallisto-2014/Sample_34967_R1.fastq/abundance.tsv \
	data/RNA-seq/kallisto-2014/Sample_34968_R1.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/67676_GTTTCG_S48_L003_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/67677_CGATGT_S49_L003_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/67678_CAGATC_S50_L004_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/67679_GTGAAA_S51_L004_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/73868_GAGTGG_S95_L006_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/73869_ATCACG_S96_L006_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/73870_TAGCTT_S97_L006_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/73871_GTGGCC_S98_L006_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/73876_ATGTCA_S4_L007_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/73884_ACTTGA_S12_L007_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/73885_CGTACG_S13_L007_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/73886_AGTTCC_S14_L007_R1_001.fastq/abundance.tsv \
	data/RNA-seq/kallisto-Run_1731/73887_GTTTCG_S15_L007_R1_001.fastq/abundance.tsv
	R -e "setwd('./src/'); source('figure_Rscripts/figure2-kallisto-post.R')"

## GO & REACTOME GSEA
./results/ECOR2HIO_24-96-RNAseq/hr24.GO-GSEA.csv ./results/ECOR2HIO_24-96-RNAseq/hr48.GO-GSEA.csv ./results/ECOR2HIO_24-96-RNAseq/hr96.GO-GSEA.csv ./results/ECOR2HIO_24-96-RNAseq/hr24.REACTOME-GSEA.csv ./results/ECOR2HIO_24-96-RNAseq/hr48.REACTOME-GSEA.csv ./results/ECOR2HIO_24-96-RNAseq/hr96.REACTOME-GSEA.csv : ./src/figure_Rscripts/figure2-GSEA.R \
	./results/ECOR2HIO_24-96-RNAseq/ECOR2_over_PBS_24hr.csv \
	./results/ECOR2HIO_24-96-RNAseq/ECOR2_over_PBS_48hr.csv \
	./results/ECOR2HIO_24-96-RNAseq/ECOR2_over_PBS_96hr.csv 
	R -e "setwd('./src/'); source('figure_Rscripts/figure2-GSEA.R')"

## Figure 2 multipanel ---------------------------------------------------------
figure2: ./figures/figure2/figure2_multipanel.pdf
./figures/figure2/figure2_multipanel.pdf : ./src/figure_Rscripts/figure2.R \
	./results/ECOR2HIO_24-96-RNAseq/complete-dataset_DESeq2-normalized-counts.csv \
	./results/ECOR2HIO_24-96-RNAseq/ECOR2_over_PBS_24hr.csv \
	./results/ECOR2HIO_24-96-RNAseq/ECOR2_over_PBS_48hr.csv \
	./results/ECOR2HIO_24-96-RNAseq/ECOR2_over_PBS_96hr.csv \
	./results/ECOR2HIO_24-96-RNAseq/hr24.GO-GSEA.csv \
	./results/ECOR2HIO_24-96-RNAseq/hr48.GO-GSEA.csv \
	./results/ECOR2HIO_24-96-RNAseq/hr96.GO-GSEA.csv \
	./results/ECOR2HIO_24-96-RNAseq/hr24.REACTOME-GSEA.csv \
	./results/ECOR2HIO_24-96-RNAseq/hr48.REACTOME-GSEA.csv \
	./results/ECOR2HIO_24-96-RNAseq/hr96.REACTOME-GSEA.csv \
	./figures/figure2/GSEA-pathways-plotted.csv \
	./data/figure1/161206_survival/survival_and_ELISA.csv
	R -e "setwd('./src/'); source('figure_Rscripts/figure2.R')"

## Figure 3 multipanel ---------------------------------------------------------
figure3: ./figures/figure3/figure3_multipanel.pdf
./figures/figure3/figure3_multipanel.pdf: ./data/figure3/dapi_ki67_edu_counts.tar.gz \
	./data/figure3/EdU_figure_cropped.png \
	./data/figure3/sox9_figure.png \
	./data/figure3/dppiv_figure.png
	tar -xzvf ./data/figure3/dapi_ki67_edu_counts.tar.gz
	R -e "setwd('./src/'); source('figure_Rscripts/figure3.R')"

## Figure 4 multipanel ---------------------------------------------------------
figure4: ./figures/figure4/figure4_multipanel.pdf
./figures/figure4/figure4_multipanel.pdf : ./src/figure_Rscripts/figure4.R \
	./data/figure4/HIO_O2_final_data.csv \
	./figures/figure4/figure4c.png
	R -e "setwd('./src/'); source('figure_Rscripts/figure4.R')"

## Generate gene counts and DE from kallisto output for Figure 5
./results/ECOR2_hypoxia_nfkb/ECOR2-HK_over_ECOR2-HK-NFKBi.csv ./results/ECOR2_hypoxia_nfkb/ECOR2-HK_over_PBS.csv ./results/ECOR2_hypoxia_nfkb/ECOR2_over_ECOR-NFKBi.csv ./results/ECOR2_hypoxia_nfkb/ECOR2_over_PBS.csv ./results/ECOR2_hypoxia_nfkb/hypoxia_over_hypoxia-NFKBi.csv ./results/ECOR2_hypoxia_nfkb/hypoxia_over_PBS.csv :
	./src/figure_Rscripts/figure5-kallisto-post.R \
	../data/RNA-seq/kallisto-Run_1731/73868_GAGTGG_S95_L006_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73869_ATCACG_S96_L006_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73870_TAGCTT_S97_L006_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73871_GTGGCC_S98_L006_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73872_AGTCAA_S99_L006_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73873_GTCCGC_S1_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73874_ACAGTG_S2_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73875_GCCAAT_S3_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73876_ATGTCA_S4_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73881_ACTGAT_S9_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73882_ATTCCT_S10_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73883_GGCTAC_S11_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73884_ACTTGA_S12_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73885_CGTACG_S13_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73886_AGTTCC_S14_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73887_GTTTCG_S15_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73888_CGATGT_S16_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73889_CAGATC_S17_L007_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73890_GTGAAA_S18_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73891_TGACCA_S19_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73892_GAGTGG_S20_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73893_ATCACG_S21_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73894_TAGCTT_S22_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73895_GTGGCC_S23_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73896_AGTCAA_S24_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73897_GTCCGC_S25_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73898_ACAGTG_S26_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73899_GCCAAT_S27_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73900_ATGTCA_S28_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73901_CTTGTA_S29_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73902_CCGTCC_S30_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73903_TTAGGC_S31_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73904_GATCAG_S32_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73905_ACTGAT_S33_L008_R1_001.fastq/abundance.h5 \
	../data/RNA-seq/kallisto-Run_1731/73906_ATTCCT_S34_L008_R1_001.fastq/abundance.h5
	R -e "setwd('./src/'); source('figure_Rscripts/figure5-kallisto-post.R')"

## Figure 5 multipanel ---------------------------------------------------------
figure5: ./figures/figure5/figure5_multipanel.pdf
./figures/figure5/figure5_multipanel.pdf ./results/ECOR2_hypoxia_nfkb/plotted-nfkb_complete-goANDreactome-results.csv : ./src/figure_Rscripts/figure5.R \
	./results/ECOR2_hypoxia_nfkb/ECOR2-HK_over_ECOR2-HK-NFKBi.csv \
	./results/ECOR2_hypoxia_nfkb/ECOR2-HK_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/ECOR2_over_ECOR-NFKBi.csv \
	./results/ECOR2_hypoxia_nfkb/ECOR2_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/hypoxia_over_hypoxia-NFKBi.csv \
	./results/ECOR2_hypoxia_nfkb/hypoxia_over_PBS.csv
	R -e "setwd('./src/'); source('figure_Rscripts/figure5.R')"

## Figure 6 multipanel ---------------------------------------------------------
figure6: ./figures/figure6/figure6_multipanel.pdf
./figures/figure6/figure6_multipanel.pdf : ./src/figure_Rscripts/figure6.R \
	./results/ECOR2_hypoxia_nfkb/ECOR2-HK_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/ECOR2_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/ECOR2-NFKBi_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/hypoxia_over_PBS.csv \
	./figures/figure6/defensin_gene_family.txt \
	./data/figure6/160518_ECOR2_BD2/160517_ECOR2_BD2_OD600.csv \
	./data/figure1/161206_survival/survival_and_ELISA.csv
	R -e "setwd('./src/'); source('figure_Rscripts/figure6.R')"

## Figure 7 multipanel ---------------------------------------------------------
figure7: ./figures/figure7/figure7_multipanel.pdf
./figures/figure7/figure7_multipanel.pdf : ./src/figure_Rscripts/figure7.R \
	./results/ECOR2_hypoxia_nfkb/ECOR2-HK_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/ECOR2_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/ECOR2-NFKBi_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/hypoxia_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/plotted-nfkb_complete-goANDreactome-results.csv \
	./figures/figure7/figure7b.png \
	./figures/figure7/figure7c.png \
	./figures/figure7/figure7d.png \
	./figures/figure7/figure7f.png
	R -e "setwd('./src/'); source('figure_Rscripts/figure7.R')"

## Figure 8 multipanel ---------------------------------------------------------
figure8: ./figures/figure8/figure8_multipanel.pdf
./figures/figure8/figure8_multipanel.pdf : ./src/figure_Rscripts/figure8.R \
	./results/ECOR2_hypoxia_nfkb/ECOR2-HK_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/ECOR2_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/ECOR2-NFKBi_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/hypoxia_over_PBS.csv \
	./results/ECOR2_hypoxia_nfkb/plotted-nfkb_complete-goANDreactome-results.csv \
	./data/figure8/results.csv \
	./data/figure8/sample_key.csv \
	./data/figure8/ECOR2_TNF_IFN_permeability.csv
	R -e "setwd('./src/'); source('figure_Rscripts/figure8.R')"
