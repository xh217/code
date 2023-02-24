maf_kataegis <- function(maf_df,min_base = min_base,distance = distance)
{

library(ClusteredMutations)

maf_df <- maf_df[,c("chr","start","ref","alt","sample")]
maf_df$chr <- gsub("chr","",maf_df$chr)

final <- c()
all_sample <- unique(maf_df$sample)
  for (i in 1:length(all_sample))
    {
       snv_sample <- maf_df[which(maf_df$sample == all_sample[i]),]
	     snv_kataegis <- features(data = snv_sample,
                       chr = chr,
                       position = start,
                       refbase = ref,
                       mutantbase = alt,
                       min = min_base,
                       max = distance)
	   if(dim(snv_kataegis)[1]!= 0)
              {
                 snv_kataegis$sample = all_sample[i]
              }
          else{
                 snv_kataegis <- data.frame() 
              }   
    final <- cbind(final,snv_kataegis)
  }

