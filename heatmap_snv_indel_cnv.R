##input df example:  sample_id  type  gene_id function_type(no space)
##                   sample_id exonic gene_id nonsynonymous

heatmap_snv_indel_cnv<-function(df)
{

df_tmp <- dcast(df,gene_id~sample_id,value.var="nonsynonymous",fun=toString)
df_tmp <- sapply(df_tmp, as.character)
df_tmp[is.na(df_tmp)] <- " "  
df_trans <- df_tmp[,-1]
rownames(df_trans) <- df_trans[,1]

library("viridis") 
function_type <- unique(df$function_type)
n_function_type <- length(function_type)
col <- viridis(10)[1:n_function_type]
names(col) <- function_type

#col = c(nonsynonymous= "brown4",stopgain="brown2", frameshift="darkblue",splicing="cyan1",CNV_gain="purple1",CNV_loss="orchid4") #customise col

alter_graphic = function(graphic = c("rect", "point"),
	                        width = 1, height = 1, 
	                        horiz_margin = unit(1, "pt"), vertical_margin = unit(1, "pt"),
	                        fill = "red", col = NA, pch = 16, ...) {

	graphic = match.arg(graphic)[1]
	if(graphic == "rect") {
		if(!is.numeric(width)) {
			stop_wrap("`width` should be nummeric.")
		}
		if(!is.numeric(height)) {
			stop_wrap("`height` should be nummeric.")
		}
		if(width != 1) {
			if(missing(horiz_margin)) {
				horiz_margin = unit(0, "pt")
			}
		}
		if(height != 1) {
			if(missing(vertical_margin)) {
				vertical_margin = unit(0, "pt")
			}
		}
		fun = function(x, y, w, h) {
			w = w*width
			h = h*height
			grid.rect(x, y, w - horiz_margin*2, h - vertical_margin*2,
				gp = gpar(fill = fill, col = col, ...))
		}
	} else if(graphic == "point") {
		fun = function(x, y, w, h) {
			grid.points(x, y, pch = pch, gp = gpar(fill = fill, col = col, ...))
		}
	}
	return(fun)
}

alter_fun.names = vector(mode = "character", length = n_function_type)
alter_fun = vector(mode = "list", length = n)

for (i in 1:n_function_type)
   {
      alter_fun[[i]] <- alter_graphic("rect", fill = col[function_type[i]])
      alter_fun.names[i] <- alter_fun.names[i]
   }   

alter_fun = list( background = alter_graphic("rect", fill = "#CCCCCC"),   
                  alter_fun)
 
heatmap_legend_param = list(title = "", at = function_type, 
                       labels = function_type )

oncoPrint(df_trans, 
          top_annotation = NULL,
          alter_fun = alter_fun,
          col = col, 
          remove_empty_rows = T,
          heatmap_legend_param = heatmap_legend_param
         )   

