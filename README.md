<p align="center">
  <img src="./misc/spaghetti-art.svg">
</p>

Spaghetti is a custom pipeline for **automatic bioinformatic analysis** of Nanopore sequencing data and **semi-automatic exploratory analysis and data visualization**. The pipeline was specifically created for the **in situ analysis of 16S rRNA gene sequences obtained by [MinION](https://nanoporetech.com/products/minion)** ([Oxford Nanopore Technologies](https://nanoporetech.com/), ONT). For that reason, Spaghetti includes tools that provide fast results and can be run on a laptop.

This pipeline was designed ad hoc for being used to characterize samples during a bioprospecting expedition to the Tabernas Desert (Almería, Spain). Read more about this in situ application here:

- Latorre-Pérez A, Gimeno-Valero H, Tanner K, Pascual J, Vilanova C and Porcar M. A round trip to the desert: in situ Nanopore sequencing informs targeted bioprospecting. *Yet to be published.*

Spaghetti consists of two different modules:

1. **Bioinformatic analysis module**: based on bash (unix). General for the analysis of 16S rRNA amplicons sequenced by ONT platforms.
2. **Exploratory analysis and data visualization module**: based on R. Specific for the characterization of the samples collected during the Tabernas Desert expedition. Some of the code can be reused for other microbial ecology analysis.

It has to be noted that **Spaghetti is not intended to be a bioinformatic tool, but a pipeline** (or a code repository) that can be reused by other users for 16S rRNA gene sequencing analysis of Nanopore data.

- [Module 1: Bioinformatic analysis](#module-1--bioinformatic-analysis)
  * [Automatic execution](#automatic-execution)
  * [References](#references)
- [Module 2: Exploratory analysis and data visualization module](#module-2-exploratory-analysis-and-data-visualization-module)
- [Why 'Spaghetti'?](#why-spaghetti)
- [Acknowledgements](#acknowledgements)
- [License](#license)
- [Citation](#citation)

# Module 1: Bioinformatic analysis

**Input**: fastq files (one fastq = one sample) generated by MinKNOW or any other basecaller.

**Output**: taxonomy and OTU-like tables that can be imported into [phyloseq](https://joey711.github.io/phyloseq/).

-----------------------
1. [Porechop](https://github.com/rrwick/Porechop) is run with default parameters for removing sequencing adapters from reads:

```{bash}
porechop -t number_of_threads -i input_file.fastq -o input_file-porechop.fastq
```

-----------------------
2. [Nanofilt](https://github.com/wdecoster/nanofilt)⁠ is used to filter reads shorter than 1,200 bp or longer than 1,800 bp.

```{bash}
cat input_file-porechop.fastq | NanoFilt -l 1200 --maxlength 1800 > input_file-porechop-nanofilt.fastq
```

-----------------------
3. Quality check is carried out with [NanoStat](https://github.com/wdecoster/nanostat)⁠.

```{bash}
NanoStat -t number_of_threads --fastq input_file-porechop-nanofilt.fastq > input_file-porechop-nanofilt-NanoStat.txt
```

-----------------------
4. Chimeras are detected and removed by using [yacrd](https://github.com/natir/yacrd) with the recommended parameters for ONT sequences

```{bash}
minimap2 -x ava-ont -g 500 -t number_of_threads input_file-porechop-nanofilt.fastq input_file-porechop-nanofilt.fastq > input_file-porechop-nanofilt.paf
yacrd -i input_file-porechop-nanofilt.paf -o input_file-porechop-nanofilt.yacrd -c 4 -n 0.4 scrubb -i input_file-porechop-nanofilt.fastq -o input_file-porechop-nanofilt.scrubb.fastq
rm *.paf
```
-----------------------
5. Filtered reads are mapped against the SILVA database (v. 138) ([Quast et al., 2013](https://academic.oup.com/nar/article/41/D1/D590/1069277))⁠, as formatted and provided by [Qiime2](https://docs.qiime2.org/2020.8/data-resources/), by using [minimap2](https://github.com/lh3/minimap2)⁠. In order to reduce minimap2’s memory usage, -K option was set to 10M, as previously suggested ([Gamaarachchi et al., 2019](https://www.nature.com/articles/s41598-019-40739-8))⁠.

The database that was used in the study can be downloaded [using this link](https://data.qiime2.org/2020.8/common/silva-138-99-seqs.qza). The path for the fasta file is: ./silva-138-99-seqs/a7432d0f-b5f7-409f-9daf-cd33db5de53f/data/dna-sequences.fasta

Taxonomy can be downloaded [using this link](https://data.qiime2.org/2020.8/common/silva-138-99-tax.qza). The path for the taxonomy file is: ./silva-138-99-tax/f12818e4-5681-434e-9af7-42b505a701d9/data/taxonomy.tsv

"./" should be replaced by the path to the folder where the files were stored.

```{bash}
# Create the database index
minimap2 -d /path/to/file/dna-sequences.mmi /path/to/file/dna-sequences.fasta

# Then map the sequences
minimap2 -x map-ont -t number_of_threads --secondary=no -K 10M /path/to/file/dna-sequences.mmi input_file-porechop-nanofilt.scrubb.fastq > input_file-porechop-nanofilt.scrubb.paf
```

The output is a [PAF file](https://github.com/lh3/miniasm/blob/master/PAF.md).

-----------------------
6. Filter the PAF files. Although secondary alignments are turned off, some query reads are actually mapped to multiple database seqs. This script filters the output file to just keep one alignment (the largest) per read. Alignments <500 pb are removed.

The script to use is [filterPAF.py](./module1/filterPAF.py).

```{bash}
filterPAF.py -i input_file-porechop-nanofilt.scrubb.paf > input_file-porechop-nanofilt.scrubb-f.paf
# Move all the PAF files (1 PAF file = 1 raw fastq file) to a new folder
mkdir filteredPAFs
mv *-f.paf filteredPAFs
```

-----------------------
7. Summarize the PAF files. This script creates a table with the number of sequences assigned to each Database's ID for each sample (assuming that 1 fastq file = 1 sample).

It's like creating an OTU table. It has to be highlighted that **Spaghetti do not work with OTUs, but with single sequences**. However, the term otu_table is used for simplicity.

The script to use is [merfePAF.py](./module1/merfePAF.py).

```{bash}
merfePAF.py -i filteredPAFs/ > otu_table.csv
```

-----------------------
8. Create a taxonomy table. This script assigns a taxonomy to every ID in the otu_table. Each ID is associated to a taxonomic category in the previously downloaded "taxonomy.tsv" file.

The script to use is [taxonomyTable.py](./module1/taxonomyTable.py).

```{bash}
taxonomyTable.py -i otu_table.csv -t /path/to/file/taxonomy.tsv > phyloseq_taxonomy.csv
```

Note that this sequence of commands will generate lots of intermediate fastq and paf files that can be removed. The most important outputs of Spaghetti are (1) the otu-like table (otu_table.csv), and (2) the taxonomy table. This tables can be directly used to import microbiome data into [phyloseq](https://joey711.github.io/phyloseq/import-data.html#phyloseq-ize_data_already_in_r)

## Automatic execution

All the commands above can be automatically run using [spaghetti.sh](./module1/spaghetti.sh). It only needs a positional argument: the folder containing all the fastq files to be analyzed. It is assumed that one fastq = one sample (= one barcode?):

```{bash}
./spaghetti.sh /path/to/fastq/files/
```

Please, note that before using [spaghetti.sh](./module1/spaghetti.sh) you should:
- Have all the tools installed and the custom scripts in your $PATH (or you can just move them to /usr/local/bin/)
- Make all the custom scripts (steps 6, 7 & 8) executable: chmod +x <fileName>
- Check the number of threads you want to use and change the "-t" parameter of the different tools
- Change the paths for the "dna-sequences.mmi" (step 5) and "taxonomy.tsv" files (step 8)
  
## References
  
This module is based on the following previous works:
  
  - [Cuscó, A., Catozzi, C., Viñes, J., Sanchez, A., & Francino, O. (2019). Microbiota profiling with long amplicons using Nanopore sequencing: full-length 16S rRNA gene and the 16S-ITS-23S of the rrn operon. F1000research, 7, 1755. doi: 10.12688/f1000research.16817.2](https://f1000research.com/articles/7-1755)
  - [Urban, L., Holzer, A., Baronas, J., Hall, M., Braeuninger-Weimer, P., & Scherm, M. et al. (2021). Freshwater monitoring by nanopore sequencing. Elife, 10. doi: 10.7554/elife.61504](https://elifesciences.org/articles/61504)
  - [Santos, A., van Aerle, R., Barrientos, L., & Martinez-Urtaza, J. (2020). Computational methods for 16S metabarcoding studies using Nanopore sequencing data. Computational And Structural Biotechnology Journal, 18, 296-305. doi: 10.1016/j.csbj.2020.01.005](https://www.sciencedirect.com/science/article/pii/S2001037019303745)

# Module 2: Exploratory analysis and data visualization module

Spaghetti was specifically designed for Latorre-Pérez *et al* (yet to be published). This work was a proof of concept for the use of in situ nanopore sequencing to characterize samples during a bioprospecting expedition. The rationale is that different sample types and/or locations can be analyzed through the journey, thus identifying which samples are more interesting depending on the objective of the bioprospecting activities. In our case, the goal was to detect those samples that showed a higher prevalence and abundance of genera that were previously described to be desiccation- and radiation-resistant. So, the data analysis was mainly -but not solely- focused on these taxa.
  
The R code for this module is provided in [spaghetti.md](./module2/spaghetti.md).
  
# Why 'Spaghetti'?

Spaghetti is a pipeline that was designed ad hoc for being used during an expedition to the Tabernas Desert (Almería, Spain). This desert was the setting chosen by [Sergio Leone](https://en.wikipedia.org/wiki/Sergio_Leone) in the 60s to shoot several of his most famous films (e.g. [Dollars Trilogy film series](https://en.wikipedia.org/wiki/Dollars_Trilogy), thus establishing the [Spaghetti Western genre](https://en.wikipedia.org/wiki/Spaghetti_Western). The name of this pipeline is a tribute to the cimenatographic legacy of this region.

<p align="center">
  <img src="./misc/tabernas-desert.jpg">
</p>

# Acknowledgements
  
I would like to thank every developer that took part in creating the tools used in the Spaghetti pipeline. I'd also like to thank the authors of the studies that inspired this pipeline ([see above](https://github.com/adlape95/Spaghetti#references)).
  
# License
  
Spaghetti is based on many open-source tools. Please check the specific licenses of each individual tool.
  
# Citation
  
- Latorre-Pérez A, Gimeno-Valero H, Tanner K, Pascual J, Vilanova C and Porcar M. A round trip to the desert: in situ Nanopore sequencing informs targeted bioprospecting. *Yet to be published.*
  
You should also cite all the tools included in the Spaghetti pipeline. Links to their GitHub pages have been included in this document. You could also find the papers cited in the "Experimental Procedures" section of our manuscript.
