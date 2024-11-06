library("QCEWAS")

ALSPAC.main.cells <- read.csv("ALSPAC.main.cells2.results.20210806.csv")

EAGeR.main.cells <- read.csv("EAGeR.main.cells.results.20210716.csv")

GenR.main.cells <- read.csv("GenR.main.cells.results.20210603.csv")

INMA.main.cells <- read.csv("INMA.main.cells.results.20210729.csv")

MoBa1.main.cells <- read.csv("MoBa1.main.cells.csv")

MoBa12.main.cells <- read.csv("MOBA2.main.cells.results.20220106.csv")

POSEIDON.main.cells <- read.csv("POSEIDON.main.cells.results.20210730.csv")

write.table(ALSPAC.main.cells, "ALSPAC.M3.txt", quote=F, row.names=F)
write.table(EAGeR.main.cells, "EAGeR.M3.txt", quote=F, row.names=F)
write.table(GenR.main.cells, "GenR.M3.txt", quote=F, row.names=F)
write.table(INMA.main.cells, "INMA.M3.txt", quote=F, row.names=F)
write.table(MoBa1.main.cells, "MoBa1.M3.txt", quote=F, row.names=F)
write.table(MoBa2.main.cells, "MoBa2.M3.txt", quote=F, row.names=F)
write.table(POSEIDON.main.cells, "POSEIDON.M3.txt", quote=F, row.names=F)

### Cohorts 
sample_list <- c("ALSPAC.M3.txt",
"EAGeR.M3.txt",
"GenR.M3.txt",
"INMA.M3.txt",
"MoBa1.M3.txt",
"MoBa2.M3.txt",
"POSEIDON.M3.txt")

## in the order you added the result files above, give the N participants
sample_N <- data.frame(file = sample_list, N = 
c(810, 354, 1162, 339, 744, 464, 293), stringsAsFactors = FALSE)


### here my microsection number has to be changed to yours (579005 --> ...)
annotation<-"/home/060089/Methylation/Annotation/450k_annotationfile_v1_2.csv"
annotation<-read.csv(annotation, header = T, strings = F)

## leave it like this 
annotation$"TARGETID"<-annotation$"IlmnID"
annotation<-annotation[,c("TARGETID", "CHR", "MAPINFO")]

## give your QC files a name (EG GENR_model1 and INMA_model1
QC_results <- EWAS_series(EWAS_files = sample_list,
output_files = c("ALSPAC_MODEL3",
"EAGeR_MODEL3",
"GENR_MODEL3",
"INMA_MODEL3",
"MoBa1_MODEL3",
"MoBa2_MODEL3",
"POSEIDON_MODEL3"),
map = annotation,
N = sample_N,
header_final_dataset = "standard",
save_final_dataset = TRUE,
high_quality_plots = TRUE,
N_plot_beta = 500000,

### if there are outliers in the result files you can filter these out or remove in this part:
threshold_outliers = c(NA, NA),
exclude_outliers = FALSE,
exclude_X = TRUE, exclude_Y = TRUE)



