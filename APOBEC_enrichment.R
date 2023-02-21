#read data
#df_mut <- read.table("/projects/TCGA_input/TCGA.ACC.somatic.maf", header = F, sep = "\t", stringsAsFactors = FALSE)
# the format of chromosomes should match that in the fasta file

tri_nuclotide_counts<-function(df_mut,bsg=BSgenome.Hsapiens.UCSC.hg38)
{
    if(exists("df_mut", mode = "list"))
      {
        mut.full <- mut.ref
      } else {
          if(file.exists(df_mut))
            {
               mut.full <- utils::read.table(df_mut, sep = "\t", header = TRUE, as.is = FALSE, check.names = FALSE)
            } else {
                stop("df_mut is neither a file nor a loaded data frame")
            }
      }	
mut=df_mut
colnames(mut) <- c("Chromosome", "Start_Position", "End_Position", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2")
mut <- mut[, c(1, 2, 4, 5, 6)]
mut <- mut[which(mut$Chromosome != "Chromosome"), ]
mut[, chr] <- factor(mut[, chr]) 
mut <- mut[which(mut[, ref] %in% c("A", "T", "C", "G") & mut[, alt] %in% c("A", "T", "C", "G")), ]
mut[, chr] <- factor(mut[, chr])
mut$context = BSgenome::getSeq(bsg, mut[, chr], as.numeric(mut[, pos]) - 1, as.numeric(mut[, pos]) + 1, as.character = T)
mut$mutcat = paste(mut[, ref], ">", mut[, alt], sep = "")
gind = grep("G", substr(mut$mutcat, 1, 1))
tind = grep("A", substr(mut$mutcat, 1, 1))
mut$std.mutcat = mut$mutcat
mut$std.mutcat[c(gind, tind)] <- gsub("G", "g", gsub("C", "c", gsub("T", "t", gsub("A", "a", mut$std.mutcat[c(gind, tind)])))) 
mut$std.mutcat[c(gind, tind)] <- gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", mut$std.mutcat[c(gind, tind)])))) 
mut$std.context = mut$context
mut$std.context[c(gind, tind)] <- gsub("G", "g", gsub("C", "c", gsub("T", "t", gsub("A", "a", mut$std.context[c(gind, tind)])))) 
mut$std.context[c(gind, tind)] <- gsub("g", "C", gsub("c", "G", gsub("t", "A", gsub("a", "T", mut$std.context[c(gind, tind)])))) 
mut$std.context[c(gind, tind)] <- sapply(strsplit(mut$std.context[c(gind, tind)], split = ""), function(str) {paste(rev(str), collapse = "")}) 
# Make the tricontext
mut$tricontext = paste(substr(mut$std.context, 1, 1), "[", mut$std.mutcat, "]", substr(mut$std.context, 3, 3), sep = "")
mut1 <- mut[, c("Tumor_Sample_Barcode", "std.mutcat", "tricontext")]
##retrive +/-20bp sequences
flank <- getSeq(bsg, mut[, chr], as.numeric(mut[, pos]) - 20, as.numeric(mut[, pos]) + 20)
flank <- as.data.frame(flank)
colnames(flank) <- "flank40"
##TC and GA counts
flank$TCcxt <- lapply(flank$flank40, function(x) {
	str_count(x, "TC")
})
flank$GAcxt <- lapply(flank$flank40, function(x) {
	str_count(x, "GA")
})
##TC[AT] and reverse context counts
flank$TCAcxt <- lapply(flank$flank40, function(x) {
	str_count(x, "TCA")
})
flank$TGAcxt <- lapply(flank$flank40, function(x) {
	str_count(x, "TGA")
})
flank$TCTcxt <- lapply(flank$flank40, function(x) {
	str_count(x, "TCT")
})
flank$AGAcxt <- lapply(flank$flank40, function(x) {
	str_count(x, "AGA")
})

##[RY]TCA and reverse context counts
flank$RTCAcxt <- lapply(flank$flank40, function(x) {
	str_count(x, "[GA]TCA")
})
flank$TGARcxt <- lapply(flank$flank40, function(x) {
	str_count(x, "TGA[GA]")
})
flank$YTCAcxt <- lapply(flank$flank40, function(x) {
	str_count(x, "[CT]TCA")
})
flank$TGAYcxt <- lapply(flank$flank40, function(x) {
	str_count(x, "TGA[CT]")
})

##G/C counts
flank$GCcxt <- lapply(flank$flank40, function(x) {
	str_count(x, "[GC]")
})
##sum corresponding reverse context
flank$TCNcxt = unlist(flank$TCcxt) + unlist(flank$GAcxt)
flank$A3cxt1 = unlist(flank$TGAcxt) + unlist(flank$TCAcxt)
flank$A3cxt2 = unlist(flank$AGAcxt) + unlist(flank$TCTcxt)
flank$YtetraCxt = unlist(flank$YTCAcxt) + unlist(flank$TGAYcxt)
flank$RtetraCxt = unlist(flank$RTCAcxt) + unlist(flank$TGARcxt)
}	
##combine and process TCW enrichment
flank1 <- flank[, c("GCcxt", "A3cxt1", "A3cxt2")]
enrich_tot <- cbind(mut1, flank1)
colnames(enrich_tot) <- c("sample", "std.mutcat", "trinuc_mut", "C_count", "TCA_count", "TCT_count")
enrich_tot$Mut_TCW <- "0"
enrich_tot$Mut_C <- "0"
enrich_tot$Con_TCW <- "0"
enrich_tot$Con_C <- "0"
# Mut_C
mutref <- data.frame(do.call("rbind", strsplit(as.character(enrich_tot$std.mutcat), ">", fixed = T)))
enrich_tot$mut_ref <- mutref[, 1]
enrich_CtoK <- enrich_tot[which(enrich_tot$std.mutcat != "C>A"), ] # Remove C>A mutations!
enrich_CtoK[which(enrich_CtoK$mut_ref == "C"), "Mut_C"] <- "1"
# Mut_TCW
enrich_CtoK[which(enrich_CtoK$trinuc_mut == "T[C>G]A"), "Mut_TCW"] <- "1"
enrich_CtoK[which(enrich_CtoK$trinuc_mut == "T[C>G]T"), "Mut_TCW"] <- "1"
enrich_CtoK[which(enrich_CtoK$trinuc_mut == "T[C>T]A"), "Mut_TCW"] <- "1"
enrich_CtoK[which(enrich_CtoK$trinuc_mut == "T[C>T]T"), "Mut_TCW"] <- "1"
# Con_C
enrich_CtoK$Con_C <- enrich_CtoK$C_count
# Con_TCW 
enrich_CtoK$Con_TCW <- enrich_CtoK$TCA_count + enrich_CtoK$TCT_count
# Aggregate and calculte enrichment score
enrich_final <- enrich_CtoK[, c("sample", "Mut_TCW", "Mut_C", "Con_TCW", "Con_C")]
enrich_final$Mut_TCW <- as.integer(enrich_final$Mut_TCW)
enrich_final$Mut_C <- as.integer(enrich_final$Mut_C)
enrich_final$Con_TCW <- as.integer(enrich_final$Con_TCW)
enrich_final$Con_C <- as.integer(enrich_final$Con_C)
enrich_final <- aggregate(enrich_final[, c(2:5)], list(enrich_final$sample), sum)
rownames(enrich_final) <- enrich_final$sample
enrich_final$sample <- NULL
enrich_final$enrich_score <- (enrich_final$Mut_TCW/enrich_final$Con_TCW)/(enrich_final$Mut_C/enrich_final$Con_C)
#Identification of samples significantly mutated by APOBEC
exe_fisher <- function(x) {
	m <- matrix(unlist(x), ncol = 2, nrow = 2, byrow = T)
	f <- fisher.test(m)
	return(as.data.frame(f$p.value))
}
enrich_matrix <- as.data.frame(enrich_final$Mut_TCW)
enrich_matrix$Mut_Denom <- enrich_final$Mut_C - enrich_final$Mut_TCW
enrich_matrix$Con_TCW <- enrich_final$Con_TCW
enrich_matrix$Con_Denom <- enrich_final$Con_C - enrich_final$Con_TCW
rownames(enrich_matrix) <- rownames(enrich_final)
colnames(enrich_matrix) <- c("Mut_TCW", "Mut_Denom", "Con_TCW", "Con_Denom")
enrich_matrix <- as.matrix(enrich_matrix)
fishers <- t(as.data.frame(apply(enrich_matrix, 1, exe_fisher)))
fishers <- as.data.frame(fishers)
enrich_final$fisher_pval <- fishers$V1
enrich_final$bh_adj_qval <- p.adjust(enrich_final$fisher_pval, method = "BH")
enrich_final$Mut_Ratio <- enrich_final$Mut_TCW/(enrich_final$Mut_C - enrich_final$Mut_TCW)
enrich_final$Con_Ratio <- enrich_final$Con_TCW/(enrich_final$Con_C - enrich_final$Con_TCW)
enrich_final[which(enrich_final$Mut_Ratio < enrich_final$Con_Ratio), "bh_adj_qval"] <- 1
enrich_final$Mut_Ratio <- NULL
enrich_final$Con_Ratio <- NULL
enrich_final$sample <- rownames(enrich_final)
rownames(enrich_final) <- NULL
return (enrich_final)
}
