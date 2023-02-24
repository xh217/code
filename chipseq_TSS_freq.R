# scoreisland_df is a dataframe with #chr1	231200	237799	205.76668758517928

chipseq_TSS_freq <- function(scoreisland_df,bsg_gene = bsg_gene,window = window...)
{
library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggplot2)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(ggprism)
library(ChIPseeker)
library(tibble)
library(tidyr)
library(ggpubr)

genome <- BSgenome.Hsapiens.UCSC.hg19
bsg_gene <- TxDb.Hsapiens.UCSC.hg19.knownGene
chrlen=seqlengths(genome)[1:23]
chrord=-paste0("chr",c(seq(1:22),"X"))
chrlen1<-rep(chrlen,each=12)

bsg_gene <- TxDb.Hsapiens.UCSC.hg19.knownGene 
bsg_gene.gr <- transcripts(txdb, columns=c("TXCHROM","TXNAME","TXSTART", "TXEND","TXSTRAND"))
bsg_gene.df=cbind.data.frame(txdb.gr$TXNAME,txdb.gr$TXCHROM, txdb.gr$TXSTART, txdb.gr$TXEND, txdb.gr$TXSTRAND)
colnames(bsg_gene)<-c("TXNAME","TXCHROM","TXSTART","TXEND","TXSTRAND")
txdb.df$TSS = ifelse( txdb.df$TXSTRAND == "+", txdb.df$TXSTART, txdb.df$TXEND )
txdb.df$TTS = ifelse( txdb.df$TXSTRAND == "-", txdb.df$TXSTART, txdb.df$TXEND)

TSS_regions <- GRanges(seqnames = Rle( txdb.df$TXCHROM ),
                      ranges = IRanges( start = txdb.df$TSS - window,
                                        end = txdb.df$TSS + window ),
                      strand = Rle( rep("*", nrow(txdb.df)) ))

TTS_regions <- GRanges(seqnames = Rle( txdb.df$TXCHROM ),
                     ranges = IRanges( start = txdb.df$TTS - window,
                                       end = txdb.df$TTS + window ),
                     strand = Rle( rep("*", nrow(txdb.df)) ))

Matrix <- read.table(scoreisland_df,header=F,sep="\t",stringsAsFactors=FALSE)
Matrix <- Matrix[,c("V1","V2","V3","V4")]
Matrix$sample <- scoreisland_df
colnames(Matrix) <- c("chr","start","end","signal","sample")
Matrix <- Matrix[which(Matrix$chr %in% chrord),]
Matrix$chr <- factor(Matrix$chr,levels=chrord)

Matrix_names = unique(Matrix$sample)       

n = length(Matrix_names)
set.names = vector(mode = "character", length = n)
set = vector(mode = "list", length = n)
  for (i in 1:n) 
    {
      Matrix_list <-  Matrix[which(Matrix$sample == Matrix_names[i]),c(1:4)]
      set[[i]] <- makeGRangesFromDataFrame(Matrix_list, keep.extra.columns=TRUE)
      set.names[i] <- Matrix_names[i]
    }
names(set) <- set.names  

tagMatrixlist_eTSS<-lapply(set, getTagMatrix, windows=TSS_regions)

plotAvgProf(tagMatrixlist_eTSS, xlim=c(-window, window), main="Repliseq peaks around TSS")+
  geom_line(aes(y = value),size=2,linetype="dashed")+
  theme_prism(base_fontface = "bold", 
  base_line_size = 1.3,
  base_size = 30) 

}
