suppressPackageStartupMessages({
    library(symphony)
    library(singlecellmethods)
    library(tidyverse)
    library(data.table)
    library(matrixStats)
    library(Matrix)
    library(plyr)
    library(dplyr)
    
    # Plotting
    library(ggplot2)
    library(ggthemes)
    library(ggrastr)
    library(RColorBrewer)
    library(patchwork)
})

fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
}