#-------------------------------------#
#1.loading packages (加载需要的包)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

#-------------------------------------#
#2 read peak file in bed format(读取bed格式的peak文件)
#2.1 bed format
#5 col: chromosome, start position, end position, peak name, score
peak <- readPeakFile('F:\\project\\OA_RA\\ATAC_part\\peak_results\\1906218RA_peaks.narrowPeak.bed')


#-------------------------------------#
#3 annotation and figure plot
#3.1 ChIP peaks coverage plot
covplot(peak, weightCol = 'V5') #all chromosome will be plot in figure
covplot(peak, weightCol = 'V5', chrs = c('chr8','chr16'), xlim = c(4.5e7, 5e7)) #specific chromosome will be plot in figure


#3.2 Profile of ChIP peaks binding to TSS regions
#3.2.1 plot figure through tagMatrix file
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000) #prepare the promoter region file
tagMatrix <- getTagMatrix(peak, windows=promoter) #make tag matrix
tagHeatmap(tagMatrix, xlim = c(-3000, 3000), color = 'red') #plot figure
#3.2.2 One step heatmap
peakHeatmap(peak, TxDb=txdb, upstream = 3000, downstream = 3000, color = 'red')


#3.3 Average Profile of ChIP peaks binding to TSS region
#3.3.1 plot figure through tagMatrix file
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000) #prepare the promoter region file
tagMatrix <- getTagMatrix(peak, windows=promoter) #make tag matrix
plotAvgProf(tagMatrix, xlim = c(-3000, 3000))
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000) # figure with confidence interval

#3.2.2 One step without tagMatrix file
plotAvgProf2(peak, TxDb = txbd, upstream = 3000, downstream = 3000,
             xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency")



#4 Peak annotation
#4.1 annotation function
peakAnno <- annotatePeak(peak, tssRegion = c(-2000,2000),
                         TxDb = txdb, annoDb = 'org.Hs.eg.db') #annoDb is optional parameter

peakAnno <- annotatePeak(peak, tssRegion = c(-2000,2000),
                         TxDb = txdb)
#4.2 figure plot function
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie = TRUE)
plotDistToTSS(peakAnno,
              title = 'Distribution of open chromosome loci\nrelative to TSs') # Visualize distribution of TF-binding loci relative to TSS


#-------------------------------------#
#5 Functional enrichment analysis
#5.1 ReactomePA for reactome pathway enrichment
library(ReactomePA)



#-------------------------------------#
#6 ChIP peak data set comparison
#6.1 Put peak files in a list
peak_1907607RA <- 'F:\\project\\OA_RA\\ATAC_part\\peak_results\\1907607RA_peaks.narrowPeak.bed'
peak_1906218RA <- 'F:\\project\\OA_RA\\ATAC_part\\peak_results\\1906218RA_peaks.narrowPeak.bed'
peak_OA2 <- 'F:\\project\\OA_RA\\ATAC_part\\peak_results\\OA2_peaks.narrowPeak'
peak_OA5 <- 'F:\\project\\OA_RA\\ATAC_part\\peak_results\\OA5_peaks.narrowPeak'
peak_RA8 <- 'F:\\project\\OA_RA\\ATAC_part\\peak_results\\RA8_peaks.narrowPeak'
peak_files <- list(peak_1907607RA=peak_1907607RA, peak_1906218RA=peak_1906218RA, 
                   peak_OA2=peak_OA2, peak_OA5=peak_OA5, peak_RA8=peak_RA8)
peak_files
#6.2 Average profiles
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000) #prepare the promoter region file
tagMatriList <- lapply(peak_files, getTagMatrix, windows = promoter)
plotAvgProf(tagMatriList, xlim = c(-3000,3000))
plotAvgProf(tagMatriList, xlim = c(-3000,3000), conf = 0.95, resample = 500, facet = 'row')
tagHeatmap(tagMatriList, xlim = c(-3000, 3000), color = NULL)
