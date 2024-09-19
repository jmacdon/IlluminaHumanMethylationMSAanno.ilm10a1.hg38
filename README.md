# IlluminaHumanMethylationMSAanno.ilm10a1.hg38

This package provides annotations for the Illumina Methylation MSA array. The data
are based on the [manifest file (csv format)](https://support.illumina.com/downloads/infinium-methylation-screening-manifest-files.html).

Presently you can install from this GitHub repo, but hopefully soon it will be part of Bioconductor.

```
library(BiocManager)
BiocManger::install("jmacdon/IlluminaHumanMethylationMSAanno.ilm10a1.hg38")

```

## Run

Currently the `minfi` package does not recognize the MSA array, so you have to 
set the annotations by hand.

```
#read intensities
MSA_example = read.metharray.exp(...)

#Annotation is assigned by hand presently
annotation(MSA_example) <- c(array = "IlluminaHumanMethylationMSA",
	                     annotation = "ilm10a1.hg38")

```
