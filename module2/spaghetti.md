---
title: "Spaghetti: In situ analysis of 16S rRNA nanopore sequencing data"
author: "Adriel Latorre-Pérez"
affiliation: "Darwin Bioprospecting Excellence"
date: "November, 2020"
version: "1.0"
output: md_document
---

```{r LibraryLoading, echo=FALSE, warning=FALSE, results=FALSE, message=FALSE}
# Load the requiered packages
library(phyloseq)
library(ggplot2)
library(plotly)
library("iNEXT")
library(ggrepel)
library(tidyr)
library(dplyr)
library("gridExtra")
```

```{r phyloseqToDF, echo=FALSE, warning=FALSE, results=FALSE, message=FALSE}
phyloseq_to_df <- function(physeq, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance"){

  # require(phyloseq)

  ## Data validation
  if(any(addtax == TRUE || sorting == "taxonomy")){
    if(is.null(phyloseq::tax_table(physeq, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }
  }

  ## Prepare data frame
  if(taxa_are_rows(physeq) == TRUE){
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), phyloseq::otu_table(physeq), stringsAsFactors = F)
  } else {
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), t(phyloseq::otu_table(physeq)), stringsAsFactors = F)
  }

  ## Check if the sample names were silently corrected in the data.frame
  if(any(!phyloseq::sample_names(physeq) %in% colnames(res)[-1])){
    if(addtax == FALSE){
      warning("Warning: Sample names were converted to the syntactically valid column names in data.frame. See 'make.names'.\n")
    }

    if(addtax == TRUE){
      stop("Error: Sample names in 'physeq' could not be automatically converted to the syntactically valid column names in data.frame (see 'make.names'). Consider renaming with 'sample_names'.\n")
    }
  }

  ## Add taxonomy
  if(addtax == TRUE){

    ## Extract taxonomy table
    taxx <- as.data.frame(phyloseq::tax_table(physeq), stringsAsFactors = F)

    ## Reorder taxonomy table
    taxx <- taxx[match(x = res$OTU, table = rownames(taxx)), ]

    ## Add taxonomy table to the data
    res <- cbind(res, taxx)

    ## Add max tax rank column
    if(addmaxrank == TRUE){

      ## Determine the lowest level of taxonomic classification
      res$LowestTaxRank <- get_max_taxonomic_rank(taxx, return_rank_only = TRUE)

      ## Reorder columns (OTU name - Taxonomy - Max Rank - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), "LowestTaxRank", phyloseq::sample_names(physeq))]

    } else {
      ## Reorder columns (OTU name - Taxonomy - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), phyloseq::sample_names(physeq))]

    } # end of addmaxrank
  }   # end of addtax

  ## Reorder OTUs
  if(!is.null(sorting)){

    ## Sort by OTU abundance
    if(sorting == "abundance"){
      otus <- res[, which(colnames(res) %in% phyloseq::sample_names(physeq))]
      res <- res[order(rowSums(otus, na.rm = T), decreasing = T), ]
    }

    ## Sort by OTU taxonomy
    if(sorting == "taxonomy"){
      taxtbl <- as.data.frame( phyloseq::tax_table(physeq), stringsAsFactors = F )

      ## Reorder by all columns
      taxtbl <- taxtbl[do.call(order, taxtbl), ]
      # taxtbl <- data.table::setorderv(taxtbl, cols = colnames(taxtbl), na.last = T)
      res <- res[match(x = rownames(taxtbl), table = res$OTU), ]
    }
  }

  ## Add OTU total abundance
  if(addtot == TRUE){
    res$Total <- rowSums(res[, which(colnames(res) %in% phyloseq::sample_names(physeq))])
  }

  rownames(res) <- NULL
  return(res)
}
```

```{r OTUtableLoading, echo=FALSE, warning=FALSE, results=FALSE}
# Load and format the (pseudo)OTU table
otufile = read.csv("./otu_table.csv", header = TRUE, sep=',')
otus = otufile[,1]
otufile = otufile[,2:length(colnames(otufile))]
rownames(otufile) = otus
otu_matrix = as.matrix(otufile)
```


```{r TaxTableLoading, echo=FALSE, warning=FALSE, results=FALSE}
# Load the taxonomy table
taxfile = read.csv("./phyloseq_taxonomy.csv", header = TRUE, sep=",")
otus = taxfile[,1]
taxfile = taxfile[,2:length(colnames(taxfile))]
rownames(taxfile) = otus
tax_matrix = as.matrix(taxfile)
```

```{r MetadataLoading, echo=FALSE, warning=FALSE, results=FALSE}
# Load the metadata
mapfile = read.csv("./Metadata.csv", header = TRUE, sep=",")
samples = mapfile[,1]
mapfile$SampleName = samples
mapfile = mapfile[,2:length(colnames(mapfile))]
rownames(mapfile) = samples
# Match the order of the sample_names
ord = match(colnames(otu_matrix), rownames(mapfile))
mapfile = mapfile[ord,]
sampledata = sample_data(mapfile)
```

```{r PhyloseqObjectCreation, echo=FALSE, warning=FALSE, results=FALSE}
# Load everything into phyloseq
OTU = otu_table(otu_matrix, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
physeq = phyloseq(OTU, TAX, sampledata)
```

# Library Size Summary

Library size by sample:

```{r SampleSums, echo=FALSE, warning=FALSE}
sample_sums(physeq)
```

```{r SampleSumsFigure, echo=FALSE, warning=FALSE}
# The same but in a Figure
data = data.frame(sample_sums(physeq), 
                  names(sample_sums(physeq)),
                  mapfile$SampleType)
colnames(data) = c("Sequences", "Sample", "Type")

fig <- plot_ly(data, x=~Sample, y = ~Sequences, color = ~Type,
  type = "bar"
)

fig
htmlwidgets::saveWidget(fig, "seqsPerSample.html")
```

Summarized library size:

```{r SAmpleSumLibrary, echo=FALSE, warning=FALSE}
summary(sample_sums(physeq))
```

-------------------------------------------------------------------------------------

```{r SpeciesCollapse, echo=FALSE, warning=FALSE}
# Collapse the phyloseq object to the species level
# -------------------------------------------------
# Explanaition: in this pipeline, 16S rRNA ONT sequences are directly mapped to a database.
# Due to the intrinsic error of ONT sequences, there is no other "easy" and "fast" way to
# set a 16S pipeline. The "problem" here is that reads that probably come from the same
# original DNA sequence will have ~94% of error, so OTU clustering it's not possible without
# transformations. This artificial discrepancies will cause that two reads coming from the
# same DNA sample will hit two different OTUs in the database. These OTUs should be taxonomically
# close (hopefully they are members of the same species). So, the solution is to collapse the
# final assignments directly into species level. This would reduce the artificial heterogeneity and
# we will be able to obtain more accurate rarefaction curves and diversity indexes.
# -------------------------------------------------

# Collapse:
physeq1 <- tax_glom(physeq, taxrank = rank_names(physeq)[7], NArm = FALSE)
# Convert into relative
physeq1_rel  = transform_sample_counts(physeq1, function(x) x / sum(x)*100 )
```

-------------------------------------------------------------------------------------

# alpha-diversity

Species level (not normalized):

```{r SpeciesAlpha, echo=FALSE, warning=FALSE, message=FALSE}
p = plot_richness(physeq1, x="SampleName", color="SampleType", measures=c("Observed","Simpson", "Shannon"), nrow = 3, sortby = "Observed") + xlab("") + ylab("") + geom_point(alpha=0.7, size = 3) + theme_light() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
p

ggsave(file = "./Figures/alfa-diversity-noNorm.svg", plot = p, device = "svg", scale = 1.5)
```

Species level (rarefied to the lower sample):

```{r SpeciesAlphaRarefied, echo=FALSE, warning=FALSE, message=FALSE}
p = plot_richness(rarefy_even_depth(physeq1, rngseed = 711), x="SampleName", color="SampleType", measures=c("Observed","Simpson", "Shannon"), nrow = 3, sortby = "Observed") + xlab("") + ylab("") + geom_point(alpha=0.7, size = 3) + theme_light() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
p

ggsave(file = "./Figures/alfa-diversity-rarefied.svg", plot = p, device = "svg", scale = 1.5)
```

Genus level (not normalized):

```{r GenusCollapse, echo=FALSE, warning=FALSE, message=FALSE}
# Collapse to Genus level:
physeq_R6 <- tax_glom(physeq, taxrank = rank_names(physeq)[6], NArm = FALSE)
# Convert into relative
physeq_R6_rel  = transform_sample_counts(physeq_R6, function(x) x / sum(x)*100 )
```

```{r GenusAlpha, echo=FALSE, warning=FALSE, message=FALSE}
p = plot_richness(physeq_R6, x="SampleName", color="SampleType", measures=c("Observed","Simpson", "Shannon"), nrow = 3, sortby = "Observed") + xlab("") + ylab("") + geom_point(alpha=0.7, size = 3) + theme_light() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
p

ggsave(file = "./Figures/GENUS-alfa-diversity-noNorm.svg", plot = p, device = "svg", scale = 1.5)
```

Genus level (rarefied to the lower sample):

```{r GenusAlphaRarefied, echo=FALSE, warning=FALSE, message=FALSE}
p = plot_richness(rarefy_even_depth(physeq_R6, rngseed = 711), x="SampleName", color="SampleType", measures=c("Observed","Simpson", "Shannon"), nrow = 3, sortby = "Observed") + xlab("") + ylab("") + geom_point(alpha=0.7, size = 3) + theme_light() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
p

ggsave(file = "./Figures/GENUS-alfa-diversity-rarefied.svg", plot = p, device = "svg", scale = 1.5)
```

-------------------------------------------------------------------------------------

# beta-diversity

PCoA - Bray Curtis - Species level:

```{r SpeciesPCoA, echo=FALSE, warning=FALSE, message=FALSE}
library(ggplot2)
# Calculate distance matrix
brayDist <- phyloseq::distance(physeq1_rel, method="bray")
# Calculate ordination
iMDS  <- ordinate(physeq1_rel, distance=brayDist, method = "PCoA")
## Make plot
# Create plot, store as temp variable, p
p <- plot_ordination(physeq1_rel, iMDS, color ="SampleType") + theme_light()
# Costumize the plot
p = p + geom_point(aes(size=2, alpha = 0.6)) + geom_text_repel(aes(label = SampleName))
p

ggsave(file = "./Figures/SPECIES-PCoA.svg", plot = p, device = "svg")
```

PCoA - Bray Curtis - Genus level:

```{r GenusPCoA, echo=FALSE, warning=FALSE, message=FALSE}
library(ggplot2)
# Calculate distance matrix
brayDist <- phyloseq::distance(physeq_R6_rel, method="bray")
# Calculate ordination
iMDS  <- ordinate(physeq_R6_rel, distance=brayDist, method = "PCoA")
## Make plot
# Create plot, store as temp variable, p
p <- plot_ordination(physeq_R6_rel, iMDS, color ="SampleType") + theme_light()
# Costumize the plot
p = p + geom_point(aes(size=2, alpha = 0.6)) + geom_text_repel(aes(label = SampleName))
p

ggsave(file = "./Figures/GENUS-PCoA.svg", plot = p, device = "svg")
```

-------------------------------------------------------------------------------------

# Top Phyla

```{r PhylumCollapse, echo=FALSE, warning=FALSE, message=FALSE}
# Collapse to Genus level:
physeq_R2 <- tax_glom(physeq1, taxrank = rank_names(physeq1)[2], NArm = FALSE)
# Convert into relative
physeq_R2_rel  = transform_sample_counts(physeq_R2, function(x) x / sum(x)*100 )
```

```{r ampvis2Loading, echo=FALSE, warning=FALSE, message=FALSE}
# Load ampvis2 package and convert the phyloseq object into a ampvis2 object

library(ampvis2)
if(!require("devtools"))
  install.packages("devtools")
#source the phyloseq_to_ampvis2() function from the gist
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

# It's mandatory to change tax_table names to match ampvis2 names
colnames(tax_table(physeq_R2_rel)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# go converty
  ampvis2_obj <- phyloseq_to_ampvis2(physeq_R2_rel)

# Change metadata names
rownames(ampvis2_obj$metadata) = ampvis2_obj$metadata[,"SampleName"]
colnames(ampvis2_obj$abund) = ampvis2_obj$metadata[,"SampleName"]
```

```{r PhylumHeatmap, echo=FALSE, warning=FALSE, message=FALSE}
p = amp_heatmap(
      data = ampvis2_obj,
      facet_by = "SampleName",
      normalise = FALSE,
      tax_show = 20,
      tax_aggregate = "Phylum",
      plot_values_size = 3,
      min_abundance = 0.000001,
      color_vector = c("gray90",
                       "whitesmoke",
                       "darkseagreen1",
                       "mediumaquamarine"),
      round = 3
    )
p

ggsave(file = "./Figures/topPhyla.svg", plot = p, device = "svg", scale = 2)
```

# Top Genera

```{r ampvis2Loading2, echo=FALSE, warning=FALSE, message=FALSE}
# Load ampvis2 package and convert the phyloseq object into a ampvis2 object

library(ampvis2)
if(!require("devtools"))
  install.packages("devtools")
#source the phyloseq_to_ampvis2() function from the gist
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

# It's mandatory to change tax_table names to match ampvis2 names
colnames(tax_table(physeq_R6_rel)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
# go converty
  ampvis2_obj <- phyloseq_to_ampvis2(physeq_R6_rel)

# Change metadata names
rownames(ampvis2_obj$metadata) = ampvis2_obj$metadata[,"SampleName"]
colnames(ampvis2_obj$abund) = ampvis2_obj$metadata[,"SampleName"]
```

```{r GeneraHeatmap, echo=FALSE, warning=FALSE, message=FALSE}
p = amp_heatmap(
      data = ampvis2_obj,
      facet_by = "SampleName",
      normalise = FALSE,
      tax_show = 20,
      tax_aggregate = "Genus",
      tax_add = "Order",
      plot_values_size = 3,
      min_abundance = 0.000001,
      color_vector = c("gray90",
                       "whitesmoke",
                       "darkseagreen1",
                       "mediumaquamarine"),
      round = 3
    )
p

ggsave(file = "./Figures/topGenera.svg", plot = p, device = "svg", scale = 2)
```


# UV-resistant bacteria:

UV-resistant and/or pigment-producing bacteria have been extensively described in the literature. So, we can know beforehand which taxa to search in this dataset. As the resolution of ONT 16S rRNA sequencing do not reach the species level (even using the entire 16S), we will search at the genus and family level.

Taxa of interest:
--

- **Genus level**: 
  - Tanner (2020): *Hymenobacter*, *Deinococcus*, *Arthrobacter*, *Alcaligenes*, *Sphingomonas*, *Curtobacterium*, *Microbacterium*, *Kineococcus*, *Methylobacterium*, *Rubrobacter*, *Modestobacter*, *Rubellimicrobium*, *Planomicrobium*, *Bacillus*, *Rhodobacter*, *Chryseobacterium*
  - Deng *et al.* (2016): *Micrococcus*
  - Paulino-Lima *et al.* (2016): *Kocuria*, *Pontibacter*, *Geodermatophilus*, *Cellulomonas*
  - Yu *et al.* (2015): *Knoellia*, *Lysobacter*, *Nocardioides*, *Paracoccus*, *Pontibacter*, *Rufibacter* and *Microvirga*
  - Etemadifar *et al.* (2015): *Exiguobacterium*
  
- **Family level** (taxonomy as in Silva 138)
  - *Hymenobacteraceae*: *Hymenobacter*, *Pontibacter*, *Rufibacter*
  - *Deinococcaceae*: *Deinococcus*
  - *Micrococcaceae*: *Arthrobacter*, *Micrococcus*, *Kocuria*
  - *Microbacteriaceae*: *Curtobacterium*, *Microbacterium*
  - *Alcaligenaceae*: *Alcaligenes*
  - *Sphingomonadaceae*: *Sphingomonas*
  - *Kineosporiaceae*: *Kineococcus*
  - *Beijerinckiaceae*: *Methylobacterium-Methylorubrum*, *Microvirga*
  - *Rubrobacteraceae*: *Rubrobacter*
  - *Geodermatophilaceae*: *Modestobacter*, *Geodermatophilus*
  - *Rhodobacteraceae*: *Rubellimicrobium*, *Rhodobacter*, *Paracoccus*
  - *Planococcaceae*: *Planomicrobium*
  - *Bacillaceae*: *Bacillus*
  - *Weeksellaceae*: *Chryseobacterium*
  - *Cellulomonadaceae*: *Cellulomonas*
  - *Intrasporangiaceae*: *Knoellia*
  - *Xanthomonadaceae*: *Lysobacter*
  - *Nocardioidaceae*: *Nocardioides*
  - *Exiguobacteraceae*: *Exiguobacterium*
  - *Trueperaceae*: *Truepera*
  
Note: *Rufibacter* and *Microvirga* are NOT in GreenGenes 13_8, but they are included in Silva 138.

(We finally use Silva 138).

------------------------------------------------

```{r SelectionUVGenus, echo=FALSE, warning=FALSE, message=FALSE}
# Let's select only the UV-resistant bacteria listed above
# I'll use the genus object, with normalization (relative proportions),
# in order to maintain the original % of the UV-bacteria selected
UV_R6_physeq = subset_taxa(physeq_R6_rel, Genus=="g__Hymenobacter" | Genus == "g__Pontibacter" | 
                             Genus == "g__Rufibacter" | Genus == "g__Deinococcus" |
                             Genus == "g__Arthrobacter" | Genus == "g__Micrococcus" |
                             Genus == "g__Kocuria" | Genus == "g__Curtobacterium" |
                             Genus == "g__Microbacterium" | Genus == "g__Alcaligenes" |
                             Genus == "g__Sphingomonas" | Genus == "g__Kineococcus" |
                             Genus == "g__Methylobacterium-Methylorubrum" | 
                             Genus == "g__Microvirga" |
                             Genus == "g__Rubrobacter" | Genus == "g__Modestobacter" |
                             Genus == "g__Geodermatophilus" | Genus == "g__Rubellimicrobium" |
                             Genus == "g__Rhodobacter" | Genus == "g__Paracoccus" |
                             Genus == "g__Planomicrobium" | Genus == "g__Bacillus" |
                             Genus == "g__Chryseobacterium" | Genus == "g__Cellulomonas" |
                             Genus == "g__Knoellia" | Genus == "g__Lysobacter" |
                             Genus == "g__Nocardioides" | Genus == "g__Exiguobacterium" |
                             Genus == "g__Truepera")
```

```{r ampvis2Loading3, echo=FALSE, warning=FALSE, message=FALSE}
# Load ampvis2 package and convert the phyloseq object into a ampvis2 object

library(ampvis2)
if(!require("devtools"))
  install.packages("devtools")
#source the phyloseq_to_ampvis2() function from the gist
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

# It's mandatory to change tax_table names to match ampvis2 names
colnames(tax_table(UV_R6_physeq)) = c("Kingdom", "Phylum", "Class", "Order", "Family", 
                                      "Genus", "Species")
# go converty
  ampvis2_obj <- phyloseq_to_ampvis2(UV_R6_physeq)

# Change metadata names
rownames(ampvis2_obj$metadata) = ampvis2_obj$metadata[,"SampleName"]
colnames(ampvis2_obj$abund) = ampvis2_obj$metadata[,"SampleName"]
```

## Percentage of UV-resistant bacteria (Genus Level)

**Not grouped:**

```{r GeneraHeatmapUV, echo=FALSE, warning=FALSE, message=FALSE}
# tax_show = 28 for showing all the UV-resistant genera
p = amp_heatmap(
      data = ampvis2_obj,
      facet_by = "SampleName",
      normalise = FALSE,
      tax_show = 29,
      tax_aggregate = "Genus",
      plot_values_size = 3,
      min_abundance = 0.000001,
      color_vector = c("gray90",
                       "whitesmoke",
                       "lightgoldenrod1",
                       "tan1",
                       "lightcoral"),
      round = 3
    )
p

ggsave(file = "./Figures/UVresistantGenera.svg", plot = p, device = "svg", scale = 2)
```

**Grouped:**

```{r GeneraBarplot, echo=FALSE, warning=FALSE, message=FALSE}
# Prepare the data.frame by joining the otu and tax table (only genus level)
bar.data = data.frame(as(otu_table(UV_R6_physeq), "matrix"), 
                      as(tax_table(UV_R6_physeq)[,"Genus"], "character"))
colnames(bar.data)[length(colnames(bar.data))] = "Genus"

# Convert the data.frame to long format
bar.data.long <- gather(bar.data, sample, abundance,
                        1:(length(colnames(bar.data)) - 1), factor_key=TRUE)

# Calculate the sum value for adding it to the fig
# Order the object by abundance
x <- bar.data.long %>%
    group_by(sample) %>%
    summarize(total = sum(abundance))
x = as.data.frame(x)
x$total = round(x$total, 2)
x = x[order(x$total, decreasing = TRUE), ]

# And plot it
bar.data.long$sample = factor(bar.data.long$sample, as.character(x$sample))
p = ggplot(data=bar.data.long, aes(x=sample, y=abundance, fill=Genus)) +
  geom_bar(stat="identity", color = "black") + theme_light() + 
  xlab("") + ylab("Relative abundance (%)")
p = p + coord_flip()
p = p + geom_text(aes(sample, total + 1, label = total, fill = NULL), data = x)

# Make the figure interactive
fig = ggplotly(p)
fig

# Save the fig
htmlwidgets::saveWidget(fig, "grouped-UVgenera.html")
```

**Richness of UV-bacteria (Presence/Absence):**

```{r GeneraUVPresenceAbsence, echo=FALSE, warning=FALSE, message=FALSE}
# Transform the data to binary
heat.data = bar.data.long
heat.data$abundance[heat.data$abundance > 0] = 1

# Summarize the richness (sum presence/absence data) 
# x$sample will serve as vector for sorting the samples by richness
x2 <- heat.data %>%
    group_by(sample) %>%
    summarize(total = sum(abundance))
x2 = as.data.frame(x2)
x2 = x2[order(x2$total, decreasing = TRUE), ]

# Plot the heatmap
heat.data$sample = factor(heat.data$sample, as.character(x2$sample))
p1 = ggplot(heat.data, aes(x = sample, y = Genus, fill = abundance)) +
  geom_tile(color="black", alpha=0.8) + xlab("") + ylab("") + theme_minimal() +
  scale_fill_gradient(low = "whitesmoke", high = "springgreen3") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1, vjust = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
p1

ggsave(file = "./Figures/UVresistantGenera_pres-abs.svg", plot = p1, device = "svg", scale = 1.5)
```

*samples has been sorted by total richness (presence of UV-resisstant bacteria)*

**Abundance + richness of UV-bacteria (all together):**

```{r GeneraUVultimateFigure, echo=FALSE, warning=FALSE, message=FALSE}
# Create the Uv-resistant bacteria abundance plot as above,
# but order by richness, not abundance (use the x2 object for that)
bar.data.long$sample = factor(bar.data.long$sample, as.character(x2$sample))
p = ggplot(data=bar.data.long, aes(x=sample, y=abundance, fill=Genus)) +
  geom_bar(stat="identity", color = "black") + theme_light() + 
  xlab("") + ylab("Relative abundance (%)")
p = p + coord_flip()
p = p + geom_text(aes(sample, total + .75, label = total, fill = NULL), data = x)
p = p + theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank())

# Create a new figure for plotting the samples sorted by richness
x2$sample = factor(x2$sample, as.character(x2$sample))
x2$pos = rep(1, length(x2$total))
p2 = ggplot(data=x2, aes(x=pos, y=sample, color=total)) +
  geom_point(size = 20) + geom_text(aes(label=total), color = "white", size = 10)
p2 = p2 + scale_color_gradient(low = "lightgoldenrod1", high = "seagreen3")
p2 = p2 + xlab("UV Bacteria") + ylab("")
p2 = p2 + theme(legend.position = "none",
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())

# Create a figure for plotting the total richness
# First estimate the richness (rarefied to the lower sample)
rich = estimate_richness(rarefy_even_depth(physeq_R6, rngseed = 711))
# Sort by UV-resistant richness
rich$pos = rep(1, length(rich$Observed))
rich$sample = rownames(rich)
rich$sample = factor(rich$sample, as.character(x2$sample))
# And do the plot
p3 = ggplot(data=rich, aes(x=pos, y=sample, color=Observed)) +
  geom_point(size = 23, shape = 18) + geom_text(aes(label=Observed), color = "white", size = 6)
p3 = p3 + scale_color_gradient(low = "lightgoldenrod1", high = "seagreen3")
p3 = p3 + xlab("Richness") + ylab("")
p3 = p3 + theme(legend.position = "none",
        axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank())

# And join the plots
grid.arrange(p, p2, p3, ncol = 3, nrow = 1, widths=c(3.5, .75, .75))
```

*Note: Samples are sorted by the total number of UV-resistant bacteria found (UV Bacteria)*
*Note 2: Rarefaction to the lower library size has been performed in order to calculate richness*

-------------------------------------------------------------------------

**OTU-like tables have been written for Species, Genera & Phylum Level.**

```{r SpeciesTable, echo=FALSE, warning=FALSE, message=FALSE}
# Species Summary Table
all = phyloseq_to_df(physeq1_rel, sorting = NULL)
# Write the table
write.table(all, row.names = FALSE, file = "./Tables/species_table.csv", quote = FALSE, sep = '\t', dec = '.')
```

```{r GenusTable, echo=FALSE, warning=FALSE, message=FALSE}
# Genus Summary Table
all = phyloseq_to_df(physeq_R6_rel, sorting = NULL)
# Write the table
write.table(all, row.names = FALSE, file = "./Tables/genus_table.csv", quote = FALSE, sep = '\t', dec = '.')
```

```{r PhylumTable, echo=FALSE, warning=FALSE, message=FALSE}
# Phylum Summary Table
all = phyloseq_to_df(physeq_R2_rel, sorting = NULL)
# Write the table
write.table(all, row.names = FALSE, file = "./Tables/phylum_table.csv", quote = FALSE, sep = '\t', dec = '.')
```

