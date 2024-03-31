# New world vulture metagenomics
 Data analysis code for the study on the gut microbiomes of new world vultures

Repository of the data and analysis procedures for the manuscript:

**[Title to be added]**
Jasper Chao, [...], Gary Graves, Antton Alberdi

## Bioinformatic procedures

Data processing to generate annotated metagenome-assembled genomes and genome count tables was conducted using the EHI bioinformatics pipeline. Data analysis procedures source from the outputs of this pipeline.

## Analysis procedures

The raw code used for data analysis is in the **Rmd** files stored in the root directory of this repository, while the bookdown-rendered webbook is available at:

[alberdilab.github.io/vulture_metagenomics](https://alberdilab.github.io/vulture_metagenomics)

While the webbook provides a user-friendly overview of the procedures, analyses can be directly reproduced using the Rmd documents. Note that the code chunks that require heavy computation have been tuned off using 'eval=FALSE'. To re-render the webbook, you can use the following code:

```r
library(bookdown)
library(htmlwidgets)
library(webshot)

render_book(input = ".", output_format = "bookdown::gitbook", output_dir = "docs")
```