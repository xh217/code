# convert vcf to maf file (header including "chr","start","ref","alt")

my_read_csv <- function(x) # x example ("/projects/vcf/example.vcf.gz")  
{
library(VariantAnnotation)
library(dplyr)
library(tidyr)
  vcf_dir = dirname(x)
  out<-readVcf(x)
    if(dim(out) == 0)
    {
  	  df = data.frame(
  	  chr = as.character(0),
  	  start = as.numeric(0),
  	  ref = as.character(0),
  	  alt = as.character(0),
  	  sample = as.character(0)
  	  )
    }
    else
    {
    df <- data.frame(info(out))
    df <- data.frame(rownames(df))
    df <- df %>% separate(colnames(df), c("chr","start","ref","alt"), "[:_/]")
    sample <- gsub(paste0(vcf_dir,"/"),"",x)
    df$sample <- sample
    df <- as.data.frame(df)
    }
  return(df)
}
