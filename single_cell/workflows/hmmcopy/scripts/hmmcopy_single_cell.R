#!/usr/bin/env Rscript

library(data.table)


suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("HMMcopy"))
suppressPackageStartupMessages(library("plyr"))
options(error=traceback)

stack_params <- function(data, paramname) {
        data = data.frame(data)
        colnames(data) = 1:length(data) - 1
        data$state = as.numeric(row.names(data)) - 1
        data = melt(data, id.vars='state', value.name='value', variable.name='iteration')
        data$parameter = paramname
        return(data)
}

format_parameter_table <- function(samp.segmented, new.params) {
        # mus - state medians
        # lambdas - state precision (inverse variance)
        # pi - state distribution
        # loglik  - likelihood values of each EM iteration

        num_iter <- ncol(samp.segmented$mus)

        loglik = stack_params(t(samp.segmented$loglik), 'loglik')
        loglik$state = NaN

        nus = stack_params(new.params$nu, 'nus')
        nus$iteration = NaN

        df.params = rbind(
                stack_params(samp.segmented$mus, 'mus'),
                stack_params(samp.segmented$lambdas, 'lambdas'),
                stack_params(samp.segmented$pi, 'pi'),
                loglik, nus)

        return(df.params)
}

error_exit_clean <- function(samp.uncorrected, chromosomes, sample_id, out_reads, out_segs, out_params, out_metrics, multiplier, error) {

        warning(paste(error, opt$tumour_file, sep=""))

        # rename space col in reads
        samp.uncorrected <- as.data.frame(samp.uncorrected)
        colnames(samp.uncorrected)[colnames(samp.uncorrected)=="space"] <- "chr"


        # uncorrected.table <- format_read_count_table(samp.uncorrected, chromosomes)
        samp.uncorrected$cell_id <- sample_id
        samp.uncorrected$cor_gc <- NA
        samp.uncorrected$cor_map <- NA
        samp.uncorrected$ideal <- FALSE
        samp.uncorrected$valid <- FALSE
        samp.uncorrected$state <- -1
        samp.uncorrected$copy <- NA
        samp.uncorrected$multiplier <- multiplier

        colorder <- c("chr","start","end","reads","gc","map","cor_gc","copy","valid","ideal","modal_curve","modal_quantile","cor_map","multiplier","state","cell_id")
        setcolorder(samp.uncorrected, colorder)

        write.table(samp.uncorrected, file=out_reads, quote=F, sep=",", col.names=T, row.names=F)

        #write colnames to the seg file
        segs <- "chr,start,end,state,median,multiplier,cell_id\n"
        cat(segs, file=out_segs)

        params <- "state,iteration,value,parameter,cell_id\n"
        cat(params, file=out_params)

        metrics_cols <- c("multiplier","MSRSI_non_integerness","MBRSI_dispersion_non_integerness",
                          "MBRSM_dispersion","autocorrelation_hmmcopy","cv_hmmcopy","empty_bins_hmmcopy",
                          "mad_hmmcopy","mean_hmmcopy_reads_per_bin","median_hmmcopy_reads_per_bin",
                          "std_hmmcopy_reads_per_bin","total_mapped_reads_hmmcopy","total_halfiness","scaled_halfiness",
                          "mean_state_mads","mean_state_vars","mad_neutral_state","breakpoints","mean_copy",
                          "state_mode","log_likelihood","true_multiplier","cell_id")

        numcols_metrics <- length(metrics_cols)
        metrics <- data.frame(matrix(c(rep.int(NA,numcols_metrics)), ncol=numcols_metrics, nrow<-1))
        colnames(metrics) <- metrics_cols
        metrics$cell_id <- sample_id
        metrics$multiplier <- multiplier
	metrics$empty_bins_hmmcopy <- 0
	metrics$total_mapped_reads_hmmcopy <- 0
	metrics$breakpoints <- 0
	metrics$state_mode <- 0
        write.table(metrics, file=out_metrics, quote=F, sep=",", col.names=T, row.names=F)
}



run_hmmcopy <- function(cell, corrected_reads_data, param, outdir, multipliers, verbose=FALSE) {

    samp.corrected <- fread(corrected_reads_data)
    samp.corrected <- data.table(start=samp.corrected$start, end=samp.corrected$end, chr=samp.corrected$chr,
                                 reads=samp.corrected$reads, gc=samp.corrected$gc, map=samp.corrected$map,
                                 cor_gc=samp.corrected$cor_gc, copy=samp.corrected$copy, valid=samp.corrected$valid, ideal=samp.corrected$ideal,
                                 modal_curve=samp.corrected$modal_curve,modal_quantile=samp.corrected$modal_quantile, cor_map=samp.corrected$cor_map)


    VALS = as.numeric(strsplit(multipliers, ",")[[1]])

    samp.corrected <- samp.corrected[order(samp.corrected$chr, samp.corrected$start), ]

    check.samp.corrected <- samp.corrected
    check.samp.corrected$copy[!check.samp.corrected$ideal] <- NaN

    #Catch and quit if no data to fit.
    if (all(is.na(check.samp.corrected$cor_gc)) | all(is.na(check.samp.corrected$copy))){

        for (VAL in VALS) {

            modal_output = file.path(outdir, VAL, sep='/')
            dir.create(modal_output, recursive=TRUE)

            out_reads <- file.path(modal_output, "reads.csv")
            out_segs <- file.path(modal_output, "segs.csv")
            out_params <- file.path(modal_output, "params.csv")
            out_metrics <- file.path(modal_output, "metrics.csv")

            err <- "Low coverage sample results in loess regression failure, unable to correct and segment"
            error_exit_clean(check.samp.corrected, chromosomes, opt$sample_id, out_reads, out_segs, out_params, out_metrics, VAL, err)
        }

        #create auto ploidy dummy output
        modal_output = file.path(outdir, '0', sep='/')
        dir.create(modal_output, recursive=TRUE)

        out_reads <- file.path(modal_output, "reads.csv")
        out_segs <- file.path(modal_output, "segs.csv")
        out_params <- file.path(modal_output, "params.csv")
        out_metrics <- file.path(modal_output, "metrics.csv")

        err <- "Low coverage sample results in loess regression failure, unable to correct and segment"
        error_exit_clean(check.samp.corrected, chromosomes, opt$sample_id, out_reads, out_segs, out_params, out_metrics, VAL, err)

        quit()

    }


    new.params <- param

    if (nrow(samp.corrected) == 0) {
        stop("INVALID INPUT")
    }

    # Initial segmentation
    seg.best <- data.frame()
    logs <- data.frame()

    best.segmented <- list()
    best.segs <- list()
    best.metrics <- list()
    best.params <- list()


    for (VAL in VALS) {

        # ROUGH
        test.corrected <- samp.corrected
        test.corrected$multiplier <- VAL
        test.corrected$copy <- test.corrected$cor_gc * VAL
        test.corrected$copy[!test.corrected$ideal] <- NaN
        samp.segmented <- HMMsegment(test.corrected, new.params, verbose = verbose, maxiter = 200)
        test.corrected$state <- samp.segmented$state - 1
        ideal <- subset(test.corrected, ideal == TRUE)

        # TWEAK
        meds <- ddply(as.data.frame(ideal), .(state), summarise, median = median(copy, na.rm = TRUE), n = length(copy))
        meds$fix <- meds$state / meds$median
        meds <- meds[order(meds$n, decreasing = TRUE), ]
        true_multiplier <- VAL * mean(subset(meds, n > 200)$fix, na.rm = TRUE)

        test.corrected$copy <- test.corrected$cor_gc * true_multiplier
        samp.segmented <- HMMsegment(test.corrected, new.params, verbose = verbose, maxiter = 200)

        # BASED 0 STATE
        test.corrected$state <- samp.segmented$state - 1
        ideal <- subset(test.corrected, ideal == TRUE)

        modal_seg <- samp.segmented$segs
        modal_seg$multiplier <- VAL
        modal_seg$state <- as.numeric(as.character(modal_seg$state)) - 1
        test.df <- as.data.frame(test.corrected)

        stats <- ddply(modal_seg, .(multiplier), summarise,
            MSRSI_non_integerness = median(abs(median - state), na.rm = TRUE)
        )

        test.df <- as.data.frame(test.corrected)
        rleseg <- rle(paste0(test.df$chr, ":", test.corrected$state))
        test.df$median <- rep(modal_seg$median, rleseg$lengths)
        test.df$halfiness <- -log2(abs(pmin(abs(test.df$median - test.df$state), 0.499) - 0.5)) - 1
        stats2 <- ddply(subset(test.df, ideal), .(multiplier), summarise,
            MBRSI_dispersion_non_integerness = median(abs(copy - state), na.rm = TRUE),
            MBRSM_dispersion = median(abs(copy - median), na.rm = TRUE),
            autocorrelation_hmmcopy = tail(acf(cor_gc, 1, na.action = na.pass, type = "correlation", plot = FALSE)$acf, 1),
            cv_hmmcopy = sd(cor_gc, na.rm = TRUE) / mean(cor_gc, na.rm = TRUE),
            empty_bins_hmmcopy = sum(reads == 0, na.rm = TRUE),
            mad_hmmcopy = mad(cor_gc, constant = 1, na.rm = TRUE),
            mean_hmmcopy_reads_per_bin = mean(reads, na.rm = TRUE),
            median_hmmcopy_reads_per_bin = median(reads, na.rm = TRUE),
            std_hmmcopy_reads_per_bin = sd(reads, na.rm = TRUE),
            total_mapped_reads_hmmcopy = sum(reads, na.rm = TRUE),
            total_halfiness = sum(halfiness, na.rm = TRUE),
            scaled_halfiness = sum(halfiness / (state + 1), na.rm = TRUE)
        )

        stats3 <- ddply(subset(test.df, ideal), .(state, multiplier), summarise,
            state_mads = mad(cor_gc, constant = 1, na.rm = TRUE),
            state_vars = var(copy, na.rm = TRUE)
        )
        stats4 <- ddply(stats3, .(multiplier), summarise,
                mean_state_mads = mean(state_mads, na.rm = TRUE),
                mean_state_vars = mean(state_vars, na.rm = TRUE)
        )

        mstats <- merge(merge(stats, stats2), stats4)
        neumad <- subset(stats3, state == 2)$state_mads
        mstats$mad_neutral_state <- ifelse(length(neumad) == 1, neumad, NA)

        mstats$breakpoints <- nrow(modal_seg) - length(unique(modal_seg$chr))
        mstats$mean_copy <- mean(ideal$copy, na.rm = TRUE)
        mstats$state_mode <- as.numeric(names(tail(sort(table(ideal$state)), 1)))
        mstats$log_likelihood <- tail(samp.segmented$loglik, 1)
        mstats$true_multiplier <- true_multiplier
        mstats$cell_id <- cell

        # HAPLOID POISON
        ones <- ideal$state == 1
        if (sum(ones) / length(ones) > 0.7) {
            mstats$scaled_halfiness <- Inf
        }

        df.params <- format_parameter_table(samp.segmented, new.params)

        # add cellid
        df.params$cell_id <- opt$sample_id
        test.corrected$cell_id <- opt$sample_id
        modal_seg$cell_id <- opt$sample_id
        mstats$cell_id <- opt$sample_id

        # rename space col in reads
        test.corrected <- as.data.frame(test.corrected)
        colnames(test.corrected)[colnames(test.corrected)=="space"] <- "chr"

        #write
        modal_output = file.path(outdir, VAL, sep='/')
        dir.create(modal_output, recursive=TRUE)
        write.table(test.corrected, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(modal_output, "reads.csv"))
        write.table(modal_seg, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(modal_output, "segs.csv"))
        write.table(mstats, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(modal_output, "metrics.csv"))
        write.table(df.params, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(modal_output, "params.csv"))

        # SAVE
        best.segmented[[VAL]] <- test.corrected
        best.segs[[VAL]] <- modal_seg
        best.metrics[[VAL]] <- mstats
        best.params[[VAL]] <- df.params

        seg.best <- rbind.fill(seg.best, data.frame(VAL, VAL, scaledpenalty = mstats$scaled_halfiness, MSRSI_non_integerness = mstats$MSRSI_non_integerness, mean_copy = mstats$mean_copy))
    }

    auto_output = file.path(outdir, "0", sep='/')
    dir.create(auto_output, recursive=TRUE)

    seg.best$red <- FALSE
    seg.best$red[which(seg.best$scaledpenalty == min(seg.best$scaledpenalty))] <- TRUE

    pick <- subset(seg.best, red)$VAL
    if (length(pick) > 1){
        pick <- pick[1]
    }

    auto_ploidy.reads <- best.segmented[[pick]]
    write.table(auto_ploidy.reads, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(auto_output, "reads.csv"))

    auto_ploidy.segs <- best.segs[[pick]]
    write.table(auto_ploidy.segs, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(auto_output, "segs.csv"))

    auto_ploidy.metrics <- best.metrics[[pick]]
    write.table(auto_ploidy.metrics, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(auto_output, "metrics.csv"))

    auto_ploidy.params <- best.params[[pick]]
    write.table(auto_ploidy.params, sep = ",", quote = FALSE, row.names = FALSE, file = file.path(auto_output, "params.csv"))
}

get_parameters <- function(str, e, mu, lambda, nu, kappa, m,eta, gamma, S) {

    str <- as.numeric(str)
    e <- as.numeric(e)
    mu <- as.numeric(strsplit(mu, ",")[[1]])
    lambda <- as.numeric(lambda)
    nu <- as.numeric(nu)
    kappa <- as.numeric(strsplit(kappa, ",")[[1]])
    m <- as.numeric(strsplit(m, ",")[[1]])

    eta <- as.numeric(eta)
    gamma <- as.numeric(gamma)
    S <- as.numeric(S)

    param <- data.frame(strength = str, e = e, 
        mu = mu, lambda = lambda, nu = nu, 
        kappa = kappa, 
        m = m, eta = eta, gamma = gamma, 
        S = S)

    return(param)
}


#=======================================================================================================================
# Command Line Options
#=======================================================================================================================
spec = matrix(c(
                "corrected_data",  "t",    1, "character", "csv file with the corrected_data",
                "sample_id",    "sample_id",    1, "character",    "specify sample or cell id",
                "outdir",      "param",    1, "character", "path to output directory",
                "param_str",      "str",    2, "double",    "optional strength parameter",
                "param_e",      "e",    2, "double",    "optional e parameter, suggested probablity of extending a segment",
                "param_mu",     "u",    2, "character", "optional mu median parameter, comma-separated list of length num_states",
                "param_l",      "l",    2, "double",    "optional lambda parameter",
                "param_nu",      "nu",    2, "double",    "optional nu parameter",
                "param_k",      "k",    2, "character", "optional kappa distribution of states parameter, comma-separated list of length num_states, should sum to 100",
                "param_m",      "p",    2, "character", "optional m median prior parameter, comma-separated list of length num_states",
                "param_eta",      "eta",    2, "character",    "optional eta parameter",
                "param_g",      "a",    2, "double",    "optional g parameter, prior shape on lambda, which is gamma distributed",
                "param_s",      "s",    2, "double",    "optional s parameter, prior scale on lambda, which is gamma distributed",
                "param_multiplier",      "mult",    2, "character",    "multiplier, start and end",
                "help",         "h",    0, "logical",   "print usage"
        ), byrow=TRUE, ncol=5);
opt = getopt(spec)

if (!is.null(opt$help)) {
    cat(getopt(spec, usage=TRUE))
    q(status=1)
}


chromosomes <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")

param <- get_parameters(opt$param_str, opt$param_e, opt$param_mu, opt$param_l, opt$param_nu, opt$param_k, opt$param_m, opt$param_eta, opt$param_g, opt$param_s)

run_hmmcopy(opt$sample_id, opt$corrected_data, param, opt$outdir, opt$param_multiplier)




