#检查是否安装了所需R包
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.15")


Package_require <- c("dplyr","ggplot2","pheatmap","cowplot",'pheatmap','tidyr','graphics','stats','tibble')

Package_install <- setdiff(Package_require,installed.packages()[, "Package"])
BiocManager::install(Package_install)
