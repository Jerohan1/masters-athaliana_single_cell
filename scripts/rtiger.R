#RTIGER not in conda so install if needed
if (!requireNamespace("RTIGER", quietly = TRUE)) devtools::install_github("rfael0cm/RTIGER") 

library(RTIGER)
setupJulia(JULIA_HOME="/faststorage/project/athaliana_single_cell/share/julia-1.12.5/bin")
#setupJulia(JULIA_HOME="/home/jool/julia-1.12.5/bin")
sourceJulia()

#si <- snakemake@input[[1]] # to mitigate snakemake not accepting a directory as input; take expected input file and take pathname.
#samples_paths  <- list.files(dirname(si), pattern = "\\.tsv$", full.names = TRUE)
samples_paths  <- list.files("../../RTIGER/col0xdb1/samples", pattern = "\\.tsv$", full.names = TRUE)
samples_paths  <- list.files("~/GenomeDK/athaliana_single_cell/workspace/RTIGER/col0xdb1/samples", pattern = "\\.tsv$", full.names = TRUE)
output_dir <- "../../RTIGER/col0xdb1"
#output_dir <- snakemake@output[[1]]
if (!dir.exists(output_dir)) 
  dir.create(output_dir)

sampleIDs <- basename(samples_paths)

# expDesign object
expDesign = data.frame(files=c(samples_paths), name=c(sampleIDs))

chr_len <- c(25768865) # from https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_028009825.2/
names(chr_len) <- c('Chr4')

myres = RTIGER(
  expDesign = expDesign,
  outputdir = output_dir,
  seqlengths = chr_len,
  rigidity = 2,
  autotune = TRUE,
  post.processing = TRUE,
  save.results = FALSE # Plotting is broken in current version, so save object to disk to plot later
)

saveRDS(myres, file = "../../RTIGER/col0xdb1/myres_object.rds")
