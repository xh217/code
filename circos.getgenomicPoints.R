circos.getgenomicPoints<- function(
	region, 
	value, 
	numeric.column = NULL, 
	sector.index = get.cell.meta.data("sector.index"),
    track.index = 1, 
    posTransform = NULL, 
	pch = par("pch"), 
	col = par("col"), 
	cex = par("cex"), 
	bg = par("bg"),
	outdir = './tmp',
	df=NULL,
	...) {
	system(paste('mkdir -p', outdir))
	nr = nrow(region)
	if(ncol(region) > 2 && inherits(region[, 1], c("character", "factor"))) {
		region = region[, -1, drop = FALSE]
	}
	
	if(is.atomic(value) && length(value) == 1) {
		value = data.frame(value = rep(value, nr))
	}
	if(is.atomic(value) && length(value) == nr) {
		value = data.frame(value = value)
	}
	if(!is.data.frame(value)) stop_wrap("`value` should be a data frame.")
	
	args = list(...)
	if(!is.null(args$.param)) {
		.param = args$.param
		if(!is.null(.param$stack)) {
			if(.param$stack && is.null(numeric.column)) {
				if(is.null(.param$jitter)) {
					value = data.frame(hline = rep(.param$i, nr))
				} else {
					value = data.frame(hline = rep(.param$i, nr) + (runif(nr) - 0.5)*abs(.param$jitter))
				}
				numeric.column = 1
			}
		} else if(!is.null(.param$numeric.column) && is.null(numeric.column)) {
			numeric.column = .param$numeric.column
		}
	}
	
	if(is.vector(value) && !is.list(value) && length(value) == 1) {
		value = data.frame(value = rep(value, nr))
		numeric.column = 1
	} else if(is.vector(value) && !is.list(value) && length(value) == nr) {
		value = data.frame(value = value)
		numeric.column = 1
	}

	if(ncol(value) == 1) numeric.column = 1
	
	if(!is.null(posTransform)) {
		region = posTransform(region)
	}
	
	if(is.null(numeric.column)) {
		numeric.column = which(as.logical(sapply(value, is.numeric)))
		if(length(numeric.column) == 0) {
			stop_wrap("Cannot find numeric column.")
		}
	}
	
	nc = length(numeric.column)
	pch = .normalizeGraphicalParam(pch, nc, nr, "pch")
	col = .normalizeGraphicalParam(col, nc, nr, "col")
	cex = .normalizeGraphicalParam(cex, nc, nr, "cex")
	bg = .normalizeGraphicalParam(bg, nc, nr, "cex")
	
		for(i in seq_len(nc)) {
			e=circos.getpoints( (region[[1]] + region[[2]])/2, value[[ numeric.column[i] ]], 
				pch = pch[i], col = col[i], cex = cex[i], bg = bg[i],
				sector.index = sector.index, track.index = 1 ) 
			print(as.data.frame(e)	)	
		}
	}

.normalizeGraphicalParam = function(x, nc, nr, name) {

	if(nc == 1) {

		if(!(length(x) == 1 || length(x) == nr)) {

			stop_wrap("The length of `", name, "` (", length(x), ") should be equal to 1 or the number of your regions (", nr, ").")

		} else if(length(x) == 1) {

			x = rep(x, nr)

		}

	} else {

		if(!(length(x) == 1 || length(x) == nc)) {

			stop_wrap("The length of `", name, "` (", length(x), ") should be equal to 1 or the number of your data column (", nc, ").")

		} else if(length(x) == 1) {

			x = rep(x, nc)

		}

	}

	return(x)

}
