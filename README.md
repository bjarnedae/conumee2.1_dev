This is a developer version (2.1). The official package will be updated soon. 

# Quick Start Guide

The conumee package provides tools for performing **co**py-**nu**mber variation (CNV) analysis using DNA **me**thylation arrays. Although the primary purpose of these arrays is the detection of genome-wide DNA methylation levels, it can additionally be used to extract useful information about copy-number alterations, e.g. in clinical cancer samples. Extracting information about CNVs from DNA methylation arrays relies on the assumption that the combined intensity values of unmethylated and methylated signals represent the copy-number state of a specific locus. Our method involves a three-step workflow: data preparation, data analysis, and output generation. Initially, we perform tangent normalization on the intensity values from a query sample to identify a unique linear combination of copy-number neutral control samples, thereby reducing technical noise. Similar to CNV analysis from genotyping arrays, this linear combination of control samples approximates the noise profile of the query sample better than any single control sample. We then calculate the log2-ratio of observed (query sample) to fitted (combined control samples) values for each probe. To further reduce technical variability, we use an adaptable heuristic to merge neighboring probes into genomic bins. The median of the normalized summed intensity values for each bin is then calculated and corrected against a baseline. Large-scale CNVs are detected using the circular binary segmentation algorithm, while focal CNVs, such as high-level amplifications and homozygous deletions, are identified in a second step. Finally, our method includes features to generate visualizations of CNVs across the genome, specific chromosomes, and individual genes, along with text-based outputs for further processing in other tools (e.g. GISTIC).

## Installation

Please make sure to install the most recent version (v2.1) of our package:

```R
remove.packages(“conumee2.0”)
devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")
```
## 1. Data preparation

### 1.1 Download test data

In this tutorial, we will analyze two Glioblastoma samples (EPICv2). For tangent normalization, we recommend using a set of at least 16 copy-number neutral control samples, ideally generated using the same experimental pipelines, and from a related biological tissue (e.g. normal human brain tissues as a control for brain tumors samples). However, we have achieved good results even with control samples that were unrelated to the query cohort. Additionally, we recommend to preprocess query and control samples in the same way (both with minfi *or* both with SeSAMe). Here, we use three control samples from the frontal cortex (450k).

```R
library("TCGAbiolinks")
library("GEOquery")

# downloading query samples (2 Glioblastoma samples from GDC, EPICv2)

query <- GDCquery(project = "CPTAC-3",
                  data.category = "DNA Methylation",
                  data.type = "Masked Intensities",
                  platform = "Illumina Methylation Epic v2",
                  barcode = c("CPT012742000", "CPT017158001"))

GDCdownload(query)

file.copy(from = list.files(paste(getwd(), "GDCdata", sep ="/"), full.names = TRUE, recursive = TRUE), # moving IDAT files into one folder
          to = paste(getwd(), "/GDCdata/CPTAC-3/DNA_Methylation/Masked_Intensities", sep = ""))
idat_dir_q <- paste(getwd(), "/GDCdata/CPTAC-3/DNA_Methylation/Masked_Intensities", sep = "")

# downloading control samples (3 normal brain samples, 450k)

GEOs <- c("GSM2403235", "GSM2403236", "GSM2403237") 
for( i in 1:3){
  getGEOSuppFiles(GEO = GEOs[i], makeDirectory = FALSE)
}

# 1. Please make sure to unzip the downloaded IDATS from the control samples and save them in a dedicated folder. 
# 2. Please create a character string named idat_dir_c with the path to that folder.
```

### 1.2 Load data

The recommended input format for Illumina 450k and EPIC arrays are MSet objects generated from IDAT files using the `minfi` package. If you would like to load data from EPICv2 arrays, mouse arrays or multiple array generations at the same time (this example), we recommend loading and preprocessing the IDAT files with the `SeSAMe` package. Both packages provide extensive functionality to preprocess the IDAT files. 

```R
library("conumee2")
library("minfi")
library("sesame")

#EPICv2 arrays, mouse arrays or IDATs from multiple arrays generations:

sdfs.q <- openSesame(idat_dir_q, prep = "QCDPB", func = NULL)
sdfs.c <- openSesame(idat_dir_c, prep = "QCDPB", func = NULL)

#__________
# commands for 450k and EPIC arrays:

RGset <- read.metharray(basenames, verbose = TRUE)
MSet <- preprocessIllumina(RGset)
```
### 1.3 Combine intensity values

The intensity values from the 'methylated' and 'unmethylated' channels are combined using the `CNV.load` function. The input can be an Mset object (minfi), a list (SeSAMe), an RnBeads object, a data.frame or matrix object containing intensities generated by GenomeStudio or downloaded from GEO (e.g., imported using read.table, which should work for most tables without reformatting). For more details, please refer to `?CNV.load`.

```R
#EPICv2 arrays, mouse arrays or IDATs from multiple arrays generations:

combined.intensities.q <- lapply(sdfs.q, totalIntensities) #lapply command for multiple samples
load.data <- CNV.load(combined.intensities.q) 

combined.intensities.c <- lapply(sdfs.c, totalIntensities) 
load.control <- CNV.load(combined.intensities.c) 

#__________
combined.intensities <- totalIntensities(sdfs) #command for one sample (not applicable for this example)
load.data <- CNV.load(combined.intensities)

# command for 450k and EPIC arrays:

load.data <- CNV.load(MSet)
```
### 1.4 Create annotation object

To begin with the CNV analysis, an annotation object, generated by the `CNV.create_anno` function, is required. This object holds information which only needs to be generated once, irrespective of the number of samples that are processed.  

- Arguments `bin_minprobes` and `bin_minsize` define the minimum number of probes per bin and the minimum bin size (default values that were optimized for 450k data are 15 and 50000, respectively). The normalized signal intensity for each bin is determined as the median log2-ratio of all the probes it contains. The genomic binning heuristic operates independently of copy-number states, ensuring that bins are consistent across samples. For baseline correction (i.e., determining the copy-number neutral state), the original bin-level log2-ratios are adjusted by a centering factor that minimizes the median absolute deviation to the baseline.  
- Argument `array_type` defines the array generation. To analyze data from multiple array generations at the same time, please choose between `overlap.1` (450k and EPIC), `overlap.2` (EPIC and EPICv2) and `overlap.3` (450k, EPIC and EPICv2). Probe annotations are loaded from the `IlluminaHumanMethylation450kanno.ilmn12.hg19` and `IlluminaHumanMethylationEPICanno.ilm10b4.hg19 packages`. For the EPICv2 and mouse arrays, probe annotations were downloaded from the manufacturer’s website and a genomic liftover was performed if necessary. In addition, information such as chromosome sizes, centromere position and gaps in the genome assembly are collected from the UCSC Genome Browser. A customized set of probes can be defined using the argument `features`.  
- Argument `genome` lets you choose between the hg19 and the hg38 probe annotation. Please note that hg38 is only supported for EPICv2 arrays.  
- Argument `exclude_regions` defines regions to be excluded (e.g. polymorphic regions, an example is given in `data(exclude_regions)`).  
- Argument `detail_regions` defines regions to be examined in detail (e.g. statistical assessment of focal copy-number status, dedicated detail plots or text output, see below). For example, detail regions can contain known oncogenes or tumor suppressor genes. These regions should either be supplied as a path to a BED file or as a GRanges object (an example is given in `data(detail_regions)`). The start and stop coordinates indicate the regions in which probes are analyzed in detail. The plotting range of detail plots are defined in a second set of start and stop coordinates in the GRanges object values (or the 7th and 8th column in the BED file). Please see `?CNV.create_anno` for more details.  

```R
data (exclude_regions) #for hg19
data(detail_regions) #for hg19
detail_regions

anno <- CNV.create_anno(array_type = "overlap.3", exclude_regions = exclude_regions, detail_regions = detail_regions) #choosing overlap.3 for analyzing EPICv2 (query) and 450k (controls) data
anno

#__________
data(detail_regions.hg38) #for hg38
```

## 2. Data processing  

### 2.1 Segmentation of the genome  

The main CNV analysis is divided into four parts:
- *Normalization*: `CNV.fit` is used to normalize a single query sample to a set of control samples through multiple linear regression. Ideally, control samples from matched normal tissues profiled within the same experiment are used, as they likely share the same technical bias. This regression analysis produces the linear combination of control samples that best matches the query sample's intensities. The log2-ratio of probe intensities between the query sample and the combined control samples is then calculated for further analysis.
  
- *Genomic Binning*: `CNV.bin` combines probes within predefined genomic bins, which are generated using `CNV.create_anno`. Intensity values are adjusted to minimize the median absolute deviation from zero across all bins, determining the copy-number neutral state.
  
- *Analysis of detail regions*: `CNV.detail` is used to analyze predefined regions in detail. This step is optional but necessary if detailed regions are to be included in plots and text files. Detail regions are defined using `CNV.create_anno`.

- *Segmentation*: `CNV.segment` segments the genome into regions with the same copy-number state. This function wraps the `CNA`, `segment`, `segments.summary`, and `segments.p` functions from the `DNAcopy` package. Default parameters are optimized for 450k data but can be modified. For more details, see `?CNV.segment`.

```R
x <- CNV.fit(load.data[1], load.control, anno)
x <- CNV.bin(x)
x <- CNV.detail(x)
x <- CNV.segment(x)
```

### 2.2 Detection of focal CNVs

Primarily designed for assessing the methylation status of specific loci, we repurpose the arrays to infer genomic copy-number changes. Notably, the microarrays are susceptible to noise and batch effects, which we aim to minimize with our analystical pipeline. However, it is important to acknowledge that methylation-based CNV detection is less sensitive than other methods and tends to provide inferential rather than definitive identification of CNVs, particularly for focal events. Therefore, we highly recommend validating positive findings with other methods, especially in clinical settings. 

The function `CNV.focal`is optional. It detects focal alterations within the predefined `detail_regions`. Complementary to this automatic calling, dedicated plots (created with `CNV.detailplot` may be helpful to decide if a focal region is significantly altered (see below). 

- Argument `sig_cgenes` is optional. It allows to assess the copy-number state for over 700 frequently altered genes from The Cancer Gene Census. Please be aware that the high number of genes may lead to false positive results.

```R
x <- CNV.focal(x)
```
## 3. Output generation 

### 3.1 Plotting functionality

The package supports multiple types of plots: 

#### Plots for single samples

- The `CNV.genomeplot` method produces plots of the complete genome or of one or multiple chromosomes. Intensity values of each bin are plotted in colored dots. Segments are shown as blue lines. If `CNV.focal` was used, significant genes are highlighted in red. See `?CNV.genomeplot` for more details.

```R
CNV.genomeplot(x)
```

![40afa67d-f399-4ca7-8d88-beadb9b2f6ad_noid_genomeplot](https://github.com/user-attachments/assets/db576f98-8848-4dec-b6e6-8c2eaa9e597f)

- The `CNV.detailplot` methods produces plots of individual detail regions, as defined in `CNV.create_anno`. Intensity values of individual probes are plotted in colored crosses. Bins are shown as black lines. Segments overlapping the detail regions are shown in blue. `CNV.detailplot_wrap` is a wrapper function that produces a single plot of all detail regions.

```R
CNV.detailplot(x, name = "EGFR")
CNV.detailplot_wrap(x)
```
![40afa67d-f399-4ca7-8d88-beadb9b2f6ad_noid_EGFR](https://github.com/user-attachments/assets/c83acb96-51ed-4811-aa05-66899c2ce40a)
![40afa67d-f399-4ca7-8d88-beadb9b2f6ad_noid_detailplot_wrap](https://github.com/user-attachments/assets/71decfcb-16e1-4f87-ae0e-c9cf822442ea)

- `CNV.plotly` creates an interactive genomeplot with annotated genes for each bin.

#### Plots for multiple samples

Please load your full cohort of query samples and follow the pipeline until `CNV.segment`.

- The `CNV.summaryplot` method converts segments from all analyzed query samples into non-overlapping, referential segments and the type of alteration (gain, loss or balanced) are summarized and visualized as percentages. The thresholds that are used for this summarization step are in line with default parameters used in GISTIC but can be adjusted by the user.

- `CNV.heatmap` generates a copy-number heatmap for all analyzed query samples.

### 3.2 Text output

Text output is generated using the `CNV.write` method. Parameter what specifies if “probes”, “bins”, “detail” , “segments”, “gistic” (for downstream processing) or “focal” (results from `CNV.focal`) should be returned. If parameter `file` is specified, the output is written into a file, otherwise a data.frame is returned. See `?CNV.write` for more details.

```R
segments <- CNV.write(x, what = "segments")
```


## 4. Contact and citation

For bug-reports, comments and feature requests please reach out via this GitHub repository. If you work with our package, please consider citing: 

Bjarne Daenekas, Eilís Pérez, Fabio Boniolo, Sabina Stefan, Salvatore Benfatto, Martin Sill, Dominik Sturm, David T W Jones, David Capper, Marc Zapatka, Volker Hovestadt, Conumee 2.0: enhanced copy-number variation analysis from DNA methylation arrays for humans and mice, Bioinformatics, Volume 40, Issue 2, February 2024, btae029, https://doi.org/10.1093/bioinformatics/btae029  

***
We thank Maximilian Leitheiser and Philipp Jurmeister for their thorough review of our code, bug fixes and helpful feedback. 









