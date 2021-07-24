<p align="center">
  <img src="./misc/spaghetti-art.svg">
</p>

Spaghetti is a custom pipeline for **automatic bioinformatic analysis** of Nanopore sequencing data and **semi-automatic exploratory analysis and data visualization**. The pipeline was specifically created for the **in situ analysis of 16S rRNA gene sequences obtained by [MinION](https://nanoporetech.com/products/minion)** ([Oxford Nanopore Technologies](https://nanoporetech.com/), ONT). For that reason, Spaghetti includes tools that provide fast results and can be run on a laptop.

This pipeline was designed ad hoc for being used to characterize samples during an expedition to the Tabernas Desert (Almería, Spain). Read more about this application here:

- Latorre-Pérez A, Gimeno-Valero H, Tanner K, Pascual J, Vilanova C and Porcar M. A round trip to the desert: in situ Nanopore sequencing informs targeted bioprospecting. *Yet to be published.*

Spaghetti consists of two different modules:

1. **Bioinformatic analysis module**: based on bash (unix). General for the analysis of 16S rRNA amplicons sequenced by ONT platforms.
2. **Exploratory analysis and data visualization module**: based on R. Specific for the characterization of the samples collected during the Tabernas Desert expedition. Some of the code can be reused for other microbial ecology analysis.

It has to be noted that **Spaghetti is not intended to be a bioinformatic tool, but a pipeline** (or a code repository) that can be used for 16S rRNA gene sequencing analysis of Nanopore data.

# Bioinformatic analysis module

1. [Porechop](https://github.com/rrwick/Porechop) is run with default parameters for removing sequencing adapters from reads:
    
```{bash}
porechop -t number_of_threads -i input_file.fastq -o input_file-porechop.fastq
```

2. [Nanofilt](https://github.com/wdecoster/nanofilt)⁠ is used to filter reads shorter than 1,200 bp or longer than 1,800 bp.

```{bash}
cat input_file-porechop.fastq | NanoFilt -l 1200 --maxlength 1800 > input_file-porechop-nanofilt.fastq
```

3. Quality check is carried out with [NanoStat](https://github.com/wdecoster/nanostat)⁠.

```{bash}
NanoStat -t number_of_threads --fastq input_file-porechop-nanofilt.fastq > input_file-porechop-nanofilt-NanoStat.txt
```

4. Chimeras are detected and removed by using [yacrd](https://github.com/natir/yacrd) with the recommended parameters for ONT sequences.

```{bash}
minimap2 -x ava-ont -g 500 -t number_of_threads input_file-porechop-nanofilt.fastq input_file-porechop-nanofilt.fastq > input_file-porechop-nanofilt.paf
yacrd -i input_file-porechop-nanofilt.paf -o input_file-porechop-nanofilt.yacrd -c 4 -n 0.4 scrubb -i input_file-porechop-nanofilt.fastq -o input_file-porechop-nanofilt.scrubb.fastq
rm *.paf
```

5. Filtered reads are mapped against the SILVA database (v. 138) ([Quast et al., 2013](https://academic.oup.com/nar/article/41/D1/D590/1069277))⁠, as formatted and provided by [Qiime2](https://docs.qiime2.org/2020.8/data-resources/), by using [minimap2](https://github.com/lh3/minimap2)⁠. In order to reduce minimap2’s memory usage, -K option was set to 10M, as previously suggested ([Gamaarachchi et al., 2019](https://www.nature.com/articles/s41598-019-40739-8))⁠.

```{bash}
minimap2
```

6. Alignments are subsequently filtered with in-house python scripts, and taxonomy and abundance tables are obtained.

# Exploratory analysis and data visualization module

# Why 'Spaghetti'?

Spaghetti is a pipeline that was designed ad hoc for being used during an expedition to the Tabernas Desert (Almería, Spain). This desert was the setting chosen by [Sergio Leone](https://en.wikipedia.org/wiki/Sergio_Leone) in the 60s to shoot several of his most famous films (e.g. [Dollars Trilogy film series](https://en.wikipedia.org/wiki/Dollars_Trilogy), thus establishing the [Spaghetti Western genre](https://en.wikipedia.org/wiki/Spaghetti_Western). The name of this pipeline is a tribute to the cimenatographic legacy of this region.

<p align="center">
  <img src="./misc/tabernas-desert.jpg">
</p>
