#!/usr/bin/env Rscript
# author: Adi Steif

# Input:
#    --tumour_file	-t 	<required> 	char	binned tumour read file (.wig)
#    --gc_file		-g	<required>	char	binned gc content file (.wig)
#    --map_file		-m 	<required>	char	binned mappability file (.wig)
#    --map_cutoff	-c 	<required>	double	optional mappability cutoff, all bins below this cutoff removed
#    --num_states	-n  <optional>  integer	number of hidden states to model, must be equal to or greater than 6
#    --param_mu     -u  <optional>  char    mu median parameter, comma-separated list of length num_states
#    --param_m      -p  <optional>  char    m median prior parameter, comma-separated list of length num_states
#    --param_k	    -k  <optional>  char	kappa distribution of states parameter, comma-separated list of length num_states, should sum to 100
#    --param_e	    -e  <optional>  double	suggested probability of staying in (extending) a segment
#    --param_g	    -a  <optional>  double	prior shape on lambda, which is gamma distributed
#    --param_s	    -s  <optional>  double	prior scale on lambda, which is gamma distributed
#    --out_dir		-o	<required>	char	path to output directory
#    --out_basename	-b	<required>	char	basic file name to which extensions will be appended
#    --help         -h	<optional>	flag	print usage

# Output:
#    [1] File with corrected read measurments (.csv)
#    [2] File with copy number segments (.csv)
#    [3] File with details of model parameter convergence and log likelihood (.csv)
#    [4] File with posterior marginals for all positions and states (.csv)







suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("HMMcopy"))
suppressPackageStartupMessages(library("plyr"))

#=======================================================================================================================
# Command Line Options
#=======================================================================================================================
spec = matrix(c(
				"corrected_data",  "t",    1, "character", "csv file with the corrected_data",
				"map_cutoff",   "c",    2, "double",    "optional mappability cutoff, all bins below this cutoff removed",
				"num_states",   "n",    2, "integer",   "optional number of hidden states to model, must be equal to or greater than 6",
				"param_mu",     "u",    2, "character", "optional mu median parameter, comma-separated list of length num_states",
				"param_m",      "p",    2, "character", "optional m median prior parameter, comma-separated list of length num_states",
				"param_k",      "k",    2, "character", "optional kappa distribution of states parameter, comma-separated list of length num_states, should sum to 100",
				"param_e",      "e",    2, "double",    "optional e parameter, suggested probablity of extending a segment",
				"param_g",      "a",    2, "double",    "optional g parameter, prior shape on lambda, which is gamma distributed",
				"param_s",      "s",    2, "double",    "optional s parameter, prior scale on lambda, which is gamma distributed",
				"param_str",      "str",    2, "double",    "optional strength parameter",
				"param_nu",      "nu",    2, "double",    "optional nu parameter",
				"param_l",      "l",    2, "double",    "optional lambda parameter",
				"param_eta",      "eta",    2, "character",    "optional eta parameter",
				"reads_output",      "r",    1, "character", "path to output directory",
				"segs_output",      "seg",    1, "character", "path to output directory",
				"params_output",      "param",    1, "character", "path to output directory",
				"sample_id",	"sample_id",	1, "character",	"specify sample or cell id",
				"post_marginals_output",      "post",    1, "character", "path to output directory",
				"auto_ploidy",         "z",    0, "logical",   "Automatically detects and corrects for ploidy, REQUIRES DIPLOID + MODAL REGRESSION",
				"help",         "h",    0, "logical",   "print usage"
		), byrow=TRUE, ncol=5);
opt = getopt(spec)

if (!is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE))
	q(status=1)
}

#=======================================================================================================================
# Helper Functions
#=======================================================================================================================
modify_param <- function(model.params, param.name, opt.param, num.states) {
	if (!is.null(opt.param)) {
		length.param <- length(strsplit(opt.param, ",")[[1]])

		if (length.param == num.states) {
			model.params[param.name] <- as.numeric(strsplit(opt.param, ",")[[1]])

		} else {
			stop(paste(paste("Invalid length for parameter ", param.name, sep=""), ". Length must be equal to num_states.", sep=""))

		}

	}
	return(model.params)
}

list_segments <- function(segs) {
	segs$chr <- as.character(segs$chr)
	segs$chr <- as.factor(segs$chr)
	segs <- split(segs, segs$chr)
}

recompute_start_end <- function(segs_chr, bin_size) {
	segs_chr <- segs_chr[,c("start", "end", "state", "median")]
	segs_chr$start <- round(segs_chr$start / bin_size) + 1
	segs_chr$end <- round(segs_chr$end / bin_size)
	colnames(segs_chr)[4] <- "median_unscaled"
	segs_chr$state <- as.character(segs_chr$state)
	segs_chr <- data.matrix(segs_chr, rownames.force=F)
}

recompute_segment_medians <- function(segs.df, samp.corrected) {
	chr <- space(samp.corrected)

	segs_list <- list_segments(segs.df)

	segs_list <- segs_list[levels(chr)]

	segs_list_rescaled <- lapply(segs_list, recompute_start_end, bin_size=samp.corrected$ranges[1]@width)

	segs_integer_median <- HMMcopy:::processSegments(segs_list_rescaled, chr, start(samp.corrected), end(samp.corrected), samp.corrected$integer_copy_scale)
}

get_bin_integer_copy_number <- function(df.row, segs) {
	within_chr <- which((as.character(segs$chr) == as.character(df.row$chr)))
	within_start <- which((segs$start <= df.row$start))
	within_end <- which((segs$end >= df.row$end))

	row_segment <- intersect(intersect(within_chr, within_start), within_end)

	if (length(row_segment)==1) {
		df.row$integer_copy_number <- segs$integer_copy_number[row_segment]
	} else {
		df.row$integer_copy_number <- NA
		warning.message <- paste("Warning! Could not find unique segment for bin:",
				df.row$chr, df.row$start, df.row$end, collapse=" ")
		warning(warning.message)
	}
}

format_read_count_table <- function(samp.corrected, chromosomes) {
	hmmcopy.table <- as.data.frame(samp.corrected)
	colnames(hmmcopy.table)[1] <- "chr"
	hmmcopy.table$chr <- factor(hmmcopy.table$chr, levels=chromosomes, ordered=T)
	hmmcopy.table <- hmmcopy.table[order(hmmcopy.table$chr),]
}

format_parameter_table <- function(samp.segmented) {
	# mus - state medians
	# lambdas - state precision (inverse variance)
	# pi - state distribution
	# loglik  - likelihood values of each EM iteration

	num_iter <- ncol(samp.segmented$mus)

	df.params <- data.frame()

	state_params <- c("mus", "lambdas")

	for (i in 1:length(state_params)) {
		df.param <- data.frame(samp.segmented[[state_params[i]]])
		colnames(df.param)[1] <- "initial"
		colnames(df.param)[2:ncol(df.param)] <- c(1:(ncol(df.param)-1))
		colnames(df.param)[ncol(df.param)] <- "final"

		if (nrow(df.param) > 1) {
			df.param$state <- c(1:nrow(df.param))
		} else {
			df.param$state <- NA
		}

		df.param$parameter <- state_params[i]

		df.params <- rbind(df.params, df.param)
	}
	df.pi <- data.frame(matrix(nrow=length(samp.segmented$pi), ncol=ncol(df.params)))
	df.pi[,ncol(df.params)-2] <- samp.segmented$pi
	df.pi[,ncol(df.params)-1] <- c(1:nrow(df.pi))
	df.pi[,ncol(df.params)] <- "pi"
	colnames(df.pi) <- colnames(df.params)
	df.params <- rbind(df.params, df.pi)

	df.loglik <- data.frame(t(samp.segmented$loglik))
	colnames(df.loglik)[1] <- "initial"
	colnames(df.loglik)[2:ncol(df.loglik)] <- c(1:(ncol(df.loglik)-1))
	colnames(df.loglik)[ncol(df.loglik)] <- "final"
	df.loglik$state <- NA
	df.loglik$parameter <- "loglik"
	df.params <- rbind(df.params, df.loglik)

	return(df.params)

}

format_posterior_marginals_table <- function(samp.corrected, samp.segmented) {
	# rho - posterior marginals (responsibilities) for each position and state

	df.bins <- as.data.frame(samp.corrected)
	colnames(df.bins)[1] <- "chr"
	df.bins <- df.bins[,c(1:4)]

	df.rho <- data.frame(t(samp.segmented$rho))
	colnames(df.rho) <- paste("state", as.character(c(1:ncol(df.rho))), sep="")

	df.rho <- cbind(df.bins, df.rho)

	df.rho$chr <- factor(df.rho$chr, levels=chromosomes, ordered=T)
	df.rho <- df.rho[order(df.rho$chr),]

	return(df.rho)

}

error_exit_clean <- function(samp.uncorrected, chromosomes, sample_id, out_reads, out_segs, out_params, out_post_marginals, error) {

	warning(paste(error, opt$tumour_file, sep=""))

	uncorrected.table <- format_read_count_table(samp.uncorrected, chromosomes)
	uncorrected.table$cell_id <- sample_id

	uncorrected.table["cor_gc"] <- NA
	uncorrected.table["cor_map"] <- NA
	uncorrected.table["ideal"] <- NA
	uncorrected.table["valid"] <- NA
	uncorrected.table["state"] <- NA
	uncorrected.table["copy"] <- NA
	uncorrected.table["integer_copy_number"] <- NA
	uncorrected.table["integer_copy_scale"] <- NA
	write.table(format(uncorrected.table, scientific=F, trim=T), file=out_reads, quote=F, sep=",", col.names=T, row.names=F)

	#write colnames to the seg file
	segs <- c("chr","start","end","state","median","integer_median","integer_copy_number")
	cat(segs, "\n", file=out_segs, sep=",")

	params <- c("initial","1","2","3","4","5","final","state","parameter","cell_id")
	cat(params, "\n", file=out_params, sep=",")

	file.create(out_post_marginals)
}

#=======================================================================================================================
# Run HMMcopy
#=======================================================================================================================
chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

#out_reads <- paste(paste(opt$out_dir, opt$out_basename, sep="/"), ".corrected_reads.csv", sep="")
#out_segs <- paste(paste(opt$out_dir, opt$out_basename, sep="/"), ".segments.csv", sep="")
#out_params <- paste(paste(opt$out_dir, opt$out_basename, sep="/"), ".parameters.csv", sep="")
#out_post_marginals <- paste(paste(opt$out_dir, opt$out_basename, sep="/"), ".posterior_marginals.csv", sep="")

out_reads <- opt$reads_output
out_segs <- opt$segs_output
out_params <- opt$params_output
out_post_marginals <- opt$post_marginals_output


samp.corrected <- read.table(opt$corrected_data, sep=',', header=TRUE)

#if correct hmmcopy fails then all corrected cols will be NA, just skip hmmcopy in that case
if (all(is.na(samp.corrected$cor_gc)) & all(is.na(samp.corrected$copy))){
		err <- "Low coverage sample results in loess regression failure, unable to correct and segment"
	error_exit_clean(samp.corrected, chromosomes, opt$sample_id, out_reads, out_segs, out_params, out_post_marginals, err)
}
samp.corrected <- RangedData(ranges = IRanges(start=samp.corrected$start, end=samp.corrected$end), space=samp.corrected$chr,
							 reads=samp.corrected$reads, gc=samp.corrected$gc, map=samp.corrected$map,
							 cor_gc=samp.corrected$cor_gc, copy=samp.corrected$copy, valid=samp.corrected$valid, ideal=samp.corrected$ideal,
							 modal_curve=samp.corrected$modal_curve,modal_quantile=samp.corrected$modal_quantile, cor_map=samp.corrected$cor_map)


if (inherits(samp.corrected, "try-error") || length((which(samp.corrected$cor.map == Inf))) > 0) {
		err <- "Low coverage sample results in loess regression failure, unable to correct and segment"
		error_exit_clean(samp.corrected, chromosomes, opt$sample_id, out_reads, out_segs, out_params, out_post_marginals, err)

} else {

	samp.corrected$copy <- samp.corrected$cor_gc

	# if mappability cutoff given, remove bins with mappability below cutoff value
	if (!is.null(opt$map_cutoff)) {
		samp.corrected$copy[samp.corrected$map < opt$map_cutoff] <- NA
	}

	# apply segmentation parameters
	default.params <- HMMsegment(samp.corrected, getparam=T, maxiter = 200)
	new.params <- default.params


	# check whether to add additional states
	if (!is.null(opt$num_states) && opt$num_states > 6) {
		while (nrow(new.params) < opt$num_states) {
			new.params <- rbind(new.params, new.params[6,])
		}

		rownames(new.params) <- NULL

		if (is.null(opt$param_mu) || is.null(opt$param_m) || is.null(opt$param_k)) {
			warning(paste("Number of states was increased, but no new values given for one of mu, m, kappa.",
						  "Parameter settings may be nonsensical.", sep=" "))
		}

	} else if (!is.null(opt$num_states) && opt$num_states < 6) {
		stop(paste(paste("Invalid value specified for num_states parameter. Must be equal to or greater than 6, but given ", opt$num_states, sep=""), ".", sep=""))
	}

	# modify model parameters (if specified)
	new.params <- modify_param(new.params, 'mu', opt$param_mu, opt$num_states)

	new.params <- modify_param(new.params, 'm', opt$param_m, opt$num_states)

	new.params <- modify_param(new.params, 'kappa', opt$param_k, opt$num_states)

	new.params <- modify_param(new.params, 'eta', opt$param_eta, opt$num_states)


	if (!is.null(opt$param_e)) {
		new.params$e <- as.numeric(opt$param_e)
	}

	if (!is.null(opt$param_str)) {
		new.params$str <- as.numeric(opt$param_str)
	}


	if (!is.null(opt$param_g)) {
		new.params$gamma <- as.numeric(opt$param_g)
	}

	if (!is.null(opt$param_s)) {
		new.params$S <- as.numeric(opt$param_s)
	}

	if (!is.null(opt$param_l)) {
		new.params$lambda <- as.numeric(opt$param_l)
	}

	if (!is.null(opt$param_nu)) {
		new.params$nu <- as.numeric(opt$param_nu)
	}

	if (is.null(opt$auto_ploidy)) {

		# segment
		tryCatch({
			samp.segmented <- HMMsegment(samp.corrected, new.params, verbose=F, maxiter = 200)
		},
			error = function(err){
				error_exit_clean(samp.corrected, chromosomes, opt$sample_id, out_reads, out_segs, out_params, out_post_marginals, as.character(err))
				quit()
		})
		samp.corrected$state <- samp.segmented$state

		# convert to integer copy number scale
		state_2_index <- which(samp.segmented$segs$state == 2)
		state_3_index <- which(samp.segmented$segs$state == 3)
		if (length(state_3_index) > 0) {

			state_3_median <- median(samp.segmented$segs$median[state_3_index])
			samp.corrected$integer_copy_scale <- samp.corrected$copy / (state_3_median / 2)

		} else if (length(state_2_index) > 0) {

			state_2_median <- median(samp.segmented$segs$median[state_2_index])
			samp.corrected$integer_copy_scale <- samp.corrected$copy / state_2_median

		} else {

			warning(paste("Sample had no State 3 or State 2 copy number calls. Unable to produce integer profile for ", opt$tumour_file, sep=""))
			samp.corrected$integer_copy_scale <- NA

		}
	} else {

		MAXMODE <- 4

		seg.best <- data.frame(MODAL = numeric(), meandiff = numeric())
		best.corrected <- list()
		best.segmented <- list()

		for (MODAL in 1:MAXMODE) {

			test.corrected <- samp.corrected
			test.corrected$copy <- test.corrected$copy * MODAL
			tryCatch({
				samp.segmented <- HMMsegment(samp.corrected, new.params, verbose=F, maxiter = 200)
			},
				error = function(err){
					error_exit_clean(samp.corrected, chromosomes, opt$sample_id, out_reads, out_segs, out_params, out_post_marginals, as.character(err))
					quit()
			})

			# BASED 0 STATE
			test.corrected$state <- samp.segmented$state - 1 

			MODAL_STATE <- as.numeric(names(sort(table(subset(test.corrected)$state), decreasing = TRUE))[1])
			modal_median <- median(subset(test.corrected, state == MODAL_STATE)$copy, na.rm = TRUE)
			test.corrected$copy <- (test.corrected$copy / modal_median)

			# PLOIDY SCALING
			N <- 5
			medsum <- ddply(as.data.frame(test.corrected), .(state), summarise, meds = median(copy, na.rm = TRUE), sd = sd(copy, na.rm = TRUE), n = length(copy))
			medsum$perc <- medsum$n / sum(medsum$n)
			medsum <- subset(medsum, !is.na(meds) & perc >= 0.01)
			medsum <- medsum[order(medsum$sd), ]
			meds <- head(medsum$meds, N)
			lowest <- data.frame(state = 0, value = Inf)
			df <- data.frame(base = c(meds, 0))
			rownames(df) <- c(head(medsum$state, N), "delta")
			for (p in 2:5) {
				multiply <- meds * p
				delta <- sum(abs(round(multiply, digits = 0) - multiply))
				df <- cbind(df, data.frame(p = c(multiply, delta)))
				colnames(df)[ncol(df)] <- p
				if (lowest$value > delta) {
					lowest <- data.frame(ploidy = p, value = delta)
				}
			}

			# EVALUATING INTEGER FIT
			test.corrected$copy <- test.corrected$copy * lowest$ploidy
			modal_seg <- samp.segmented$segs
			modal_seg$state <- as.numeric(modal_seg$state) - 1
			modal_seg$median <- modal_seg$median * (lowest$ploidy / modal_median)
			modal_seg$closest_int <- abs(modal_seg$median - as.numeric(as.character(modal_seg$state)))
			modal_seg$halfway <- abs(modal_seg$closest_int - 0.5)
			modal_seg <- modal_seg[order(modal_seg$halfway), ]

			halfway <- subset(modal_seg, state < 4 & halfway <= 0.1) # LOL
			halfway$penalty <- abs(halfway$halfway - 0.1)
			halfway$width <- halfway$end - halfway$start
			halfway$score <- halfway$penalty
			halfiness <- sum(halfway$score)

			medians <- ddply(modal_seg, .(state), summarise, med = median(median))
			medians$diff <- abs(medians$state - medians$med)
			meandiff <- mean(subset(medians, state %in% 0:4)$diff)

			# SAVING RESULTS
			seg.best <- rbind.fill(seg.best, data.frame(MODAL, meandiff, halfiness, delta = lowest$value))
			best.corrected[[MODAL]] <- test.corrected
			best.segmented[[MODAL]] <- samp.segmented

		}

		seg.best <- seg.best[order(seg.best$meandiff, seg.best$MODAL, decreasing = FALSE), ]
		samp.corrected <- best.corrected[[seg.best$MODAL[1]]]
		samp.segmented <- best.segmented[[seg.best$MODAL[1]]]
	}

	# recompute segment medians
	segs.integer.medians <- recompute_segment_medians(samp.segmented$segs, samp.corrected)
	colnames(segs.integer.medians)[ncol(segs.integer.medians)] <- "integer_median"

	# round segment medians to integer copy number
	segs.integer.medians$integer_copy_number <- round(segs.integer.medians$integer_median)

	# format corrected read count table
	corrected.table <- format_read_count_table(samp.corrected, chromosomes)

	# add integer copy number to read count table
	corrected.table <- adply(corrected.table, 1, get_bin_integer_copy_number, segs=segs.integer.medians)
	colnames(corrected.table)[ncol(corrected.table)] <- "integer_copy_number"
		corrected.table$cell_id <- opt$sample_id
	# output read count table\

		# names(corrected.table)[names(corrected.table) == 'cor.map'] <- 'cor_map'
		# names(corrected.table)[names(corrected.table) == 'cor.gc'] <- 'cor_gc'
		#print(colnames(corrected.table))
	write.table(format(corrected.table, scientific=F, trim=T), file=out_reads, quote=F, sep=",", col.names=T, row.names=F)

	# format and output segment table
	segments.table <- merge(samp.segmented$segs, segs.integer.medians, sort=F)
	segments.table$chr <- factor(segments.table$chr, levels=chromosomes, ordered=T)
	segments.table <- segments.table[order(segments.table$chr),]
		segments.table$cell_id <- opt$sample_id
	write.table(format(segments.table, scientific=F, trim=T), file=out_segs, quote=F, sep=",", col.names=T, row.names=F)

	# format and output parameter and posterior marginal tables
	df.params <- format_parameter_table(samp.segmented)

	#add nus data to params output
	nus <- data.frame(matrix(ncol = ncol(df.params), nrow=nrow(new.params)))
	colnames(nus) <- colnames(df.params)
	nus$final <- new.params$nu
	nus$parameter <- "nus"
	nus$state <- rownames(new.params)

	df.params <- rbind(df.params, nus)
	df.params$cell_id <- opt$sample_id

	write.table(format(df.params, scientific=F, trim=T), file=out_params, quote=F, sep=",", col.names=T, row.names=F)

	df.marginals <- format_posterior_marginals_table(samp.corrected, samp.segmented)
		df.marginals$cell_id <- opt$sample_id
	write.table(format(df.marginals, scientific=F, trim=T), file=out_post_marginals, quote=F, sep=",", col.names=T, row.names=F)

}

# hTERT
#opt$tumour_file <- "/share/scratch/asteif_temp/single_cell_indexing/hmmcopy/PX0281/SA040-P17-2/2016-02-01_bin_200000_e_default/tmp/SA040-01818.tumour.wig"
#opt$out_dir <- "/share/lustre/asteif/projects/single_cell_indexing/test/hmmcopy_test"
#opt$out_basename <- "SA040-01818"


#opt$gc_file <- "/share/lustre/asteif/applications/HMMcopy/HMMcopy/data/GRCh37-lite.gc.ws_200000.wig"
#opt$map_file <- "/share/lustre/asteif/applications/HMMcopy/HMMcopy/data/GRCh37-lite.map.ws_125_to_200000.wig"
#opt$map_cutoff <- 0.9
#opt$num_states <- 7
#opt$param_mu <- "0,0.5,1.0,1.5,2.0,2.5,3.0"
#opt$param_m <- "0,0.5,1.0,1.5,2.0,2.5,3.0"
#opt$param_k <- "25,50,800,50,25,25,25"
#opt$param_e <- 0.9999999
#opt$param_g <- 3
#opt$param_s <- 0.1592395
