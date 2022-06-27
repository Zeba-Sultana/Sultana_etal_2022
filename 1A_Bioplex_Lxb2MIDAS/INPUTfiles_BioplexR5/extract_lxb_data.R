# .lxb files extraction script
#' @param plate An array describing the organisation of the plate, must be defined in the R session
#' @param lxb_folder The folder containing the lxb files
library("lxb")

# join function to create a string from a list of string
join = function(string_list, sep=" ") {
    first=TRUE
    for (element in string_list) {
        if (element != "") {
            if(first) {
                full_string = element;
                first = FALSE;
            } else {
                full_string = paste(full_string, element, sep=sep)
            }
        }
    }
    return(full_string)
}
# Personnalized boxplot
pbox = function(data, main="", sub="", xlab, ylab) {
    boxplot(data, main=main, sub=sub, xlab=xlab, ylab=ylab, range=0)
}
# Converging boostrap, with a minimum number of run to ensure some initial variation
#' @return a vector or length 3 with the mean, and the 'error' confidence interval of this mean. 5% by default (0.025 percentile on each side)
bootstrapping = function(values, func=median, error=0.05, min_run=100, max_run=1000*1000, EPS=0.001) {
    bootstrap = c()
    for (i in 1:min_run) {
        med = func(sample(values, replace=T))
        bootstrap = c(bootstrap, med)
    }

    old_mean = mean(bootstrap)
    new_mean = old_mean + EPS
    while(abs(old_mean - new_mean) >= EPS && length(bootstrap) < max_run)
    {
        #for (i in 1:min_run) {
            #med = func(sample(values, replace=T))
            #bootstrap = c(bootstrap, med)
        #}
        bootstrap = c(bootstrap, sapply(1:min_run, function(X) { func(sample(values, replace=T)) }))
        
        old_mean = new_mean
        new_mean = mean(bootstrap)
    }
    bootstrap = sort(bootstrap)

    return(c(mean(bootstrap), bootstrap[floor(error/2 * length(bootstrap))], bootstrap[ceiling((1-error/2) * length(bootstrap))]))
}

# Get the association between wells and treatments
extractExperimentalSetup <- function(plate) {
    treatments = list()
    for (row in rownames(plate)) {
        for (col in colnames(plate)) {
            well = paste0(row, col)
            if (plate[row, col] != "") {
                treatments[[well]] = plate[row, col]
            }
        }
    }
    wells_per_treatment = list()
    blanks = c(); controls = c(); externals=c()
    for (well in names(treatments)) {
        if (grepl("[bB]lank", treatments[[well]])) { # Blank wells
            blanks = c(blanks, well)
        } else if (grepl( "control", tolower(treatments[[well]]) )) {
            controls = c(controls, well)
        } else if (length(treatments[[well]]) == 0) {
            externals = c(externals, well)
        } else if (!grepl("[Kk][Oo]", treatments[[well]])) {
            wells_per_treatment[[ treatments[[well]] ]] = c(wells_per_treatment[[ treatments[[well]] ]], well)
        }
    }
    wells_per_treatment_with_blank = wells_per_treatment;
    wells_per_treatment_with_blank[["blank"]] = blanks;
# Extract stimulators and perturbators
    perturbators = unique(unlist(strsplit(unique(unlist(treatments)), "\\+")))
    perturbators = perturbators[perturbators!="control" & perturbators!="blank"]
    stimulators = perturbators[!grepl("i$", perturbators)]
    inhibitors = perturbators[grepl("i$", perturbators)]

    return(list(treatments=treatments, wptb=wells_per_treatment_with_blank, wpt=wells_per_treatment, stimulators=stimulators, inhibitors=inhibitors, controls=controls, blanks=blanks, externals=externals))
}

#' Process the bead dataset
#'
#' @param lxb_dataset The lxb dataset
#' @param region Beads IDs used
#' @param analyte Analyte corresponding to the analysis.
#' @param controls Wells containing control conditions.
#' @param blanks Wells containing no cell.
#' @param wells_per_treatment_with_blank List of wells where each treatment is applied. It is a list where each entry is named after a treatment and contains a vector of well names.
#' @param externals Wells to exclude from the processing
#' @param PLOT Whether data should be plotted
readBeads <- function(lxb_dataset, region, analyte, controls, blanks, wells_per_treatment_with_blank, plot_dir="", externals=c(""), PLOT=F) {
    wells_per_treatment = wells_per_treatment_with_blank
    wells_per_treatment[["blank"]] = NULL # Suppress this entry
    # Collection of samples caracteristics
    cell_data = c()
    cell_variation = c()
    # Local collection of the histogram data
    normality = c()
    dist_beads = c()
    invalids = c()
    #wells_variation = data.frame(row.names=names(lxb_dataset))
    wells_variation = c()
    for (bead in region) {
        antibody = analyte[which(region == bead)]
        dist_beads[[antibody]] = list()
        bead_carac = c(); # Value of the antibody in each well
        bead_blank = c(); # Value of the blank for this antibody
        bead_control = c(); # Value of the control for this antibody
        median_cv = c(); # Collect the sample CV from the bootstrap
        if(PLOT) { pdf(paste(plot_dir, gsub("/", "-", gsub(" ", "_", analyte[region==bead])), "_", basename(plot_dir), ".pdf", sep="")) }
        for (well in names(lxb_dataset)) {
            if (!(well %in% externals)) {
                # Extract the information (bootstrapped median + error)
                processing = processReadout(lxb_dataset, bead, well, PLOT)
                if (length(processing$values) <  1) { # If there are no beads, do not use the data
                    #warning(paste0("Well ", well, " has 0 beads of type ", analyte[which(region==bead)]))
                    invalids = rbind(invalids, c(well, analyte[which(region==bead)]))
                }
                # Gather the information with a convenient indexing
                dist_beads[[antibody]][[well]] = processing$values
                normality = c(normality, processing$shapiro)
                # Bead distribution, per cell line and per antibody
                bead_carac[well] = processing$bootstrap[1]
                median_cv[well] = processing$cv
            }
        }
        if (PLOT) { dev.off() }

        variations = computeVariation(wells_per_treatment_with_blank, bead_carac, dist_beads[[antibody]], median_cv, blanks, controls)
        wells_variation = cbind(wells_variation, median_cv)
        cell_variation = cbind(cell_variation, variations$cv)
        cell_data = cbind(cell_data, bead_carac)
    }
    colnames(wells_variation) = analyte
    colnames(cell_variation) = analyte
    colnames(cell_data) = analyte
    if (length(invalids) > 0) { warning(paste( "Some wells have no beads for some readouts:", paste0(unique(invalids[,1]), collapse=", ") )) }
    return(list(datas=cell_data, cv=cell_variation, controls=controls, blanks=blanks, wpt=wells_per_treatment_with_blank, invalid_wells=invalids, normality=normality, wells_variation=wells_variation))
}

#' Extract data for a bead number in a well from an lxb dataset
#' 
#' Select valid data for the specified bead of the specified well and perform a bootstrap to extract a robust mean an the variation
#' @param dataset An full lxb dataset
#' @param bead The bead number
#' @param well The well identifier
processReadout <- function(dataset, bead, well, PLOT=FALSE) {
    data = dataset[[well]]

    # Select valid measurements for the selected bead and get rid of outlayers for the plot
    # CL = Classification channel and RP = Readout value
    valids = which(data[,"Bead Valid"]==1 & data[,"Bead ID"]==bead & data[, "RP1S Valid"]==1 & data[, "RP1L Valid"]==1 & data[, "CL1 Valid"]==1 & data[, "CL2 Valid"]==1 )
    values = data[valids, "RP1 Value"]
    if (length(valids) < 3) { # 3 values are required for the shapiro test
        return(list(values=values, shapiro=c(), cv=0, bootstrap=c(0, 0, 0)))
    }
    # Bootstrap the median on all the data
    bst = bootstrapping(values, error=0.05)
    #run2 = c(run2, bootstrapping(data[valids, "RP1 Value"])[1])
    median_cv = (bst[3] - bst[2]) / (bst[1] * 2 * 1.96) # 2 * 1.96 to have one sd under normality assumption
    # Plot the histogram for each bead in the well with a normal density on top of it
    if(PLOT) {# False for quick computation
        plot_values = data[valids[data[valids, "RP1 Value"] < ( median(data[valids,"RP1 Value"] + 3*sd(data[valids, "RP1 Value"])) )], "RP1 Value"]
        par(col="black")
        hist(plot_values, 1 + round(length(plot_values)/5), main=paste(analyte[which(region==bead)], "for well", well), prob=TRUE)

        med = median(plot_values);
        dev = sd(plot_values)
        range = med + 4*dev;
        nb_points = 10 * 1000;
        par(col="red")
        lines(seq(0, range, range/nb_points), dnorm(seq(0, range, range/nb_points), mean=med, sd=dev))
        # Bootstrap results
        lines(rep(bst[1], 101), 0:100/1000, col="blue", lty=2, lwd=1.5)
        lines(rep(bst[2], 101), 0:100/1000, col="blue", lty=2, lwd=0.7)
        lines(rep(bst[3], 101), 0:100/1000, col="blue", lty=2, lwd=0.7)
        mtext(paste(length(plot_values), "beads"), 3)
    }
    # Normality test, collection per cell line and per antibody
    normal = tryCatch(shapiro.test(values)$statistic, error=function(X) {return(1)})

    return(list(values=values, shapiro=normal, cv=median_cv, bootstrap=bst))
}

#' Compute the coefficient of variation for one readout
#'
#' Compute the coefficient of variation for one readout taking into account the variation of the for each well and the variation between samples
#' @param wells_per_treatment List of vectors, each entry of the list is named after a treatment and contains a list of wells coordinates
#' @param bead_carac Values for all the wells for one bead
#' @param wilcox_sample List of the distribution of values of the bead in each well
computeVariation <- function (wells_per_treatment, bead_carac, wilcox_sample, median_cv, blanks, controls, default_cv=0.3, min_cv=0.1) {
    # Get CVs from replicated conditions that are not too close to blank (i.e where CV != noise) to have a global noise level for the readout
    replicates_cv = c()
    for (wells in wells_per_treatment) {
        collect = c()
        for (well in wells) {
            # Without blank, we don't know the background noise level so we collect all
            if (length(blanks) == 0 || ( !is.nan(mean(bead_carac[blanks], na.rm=TRUE)) && bead_carac[well] > 2 * mean(bead_carac[blanks], na.rm=TRUE) )) {
                collect = c(collect, well)
            }
        }
        if (length(collect) > 1) {
            replicates_cv = c( replicates_cv, sd(bead_carac[wells]) / mean(bead_carac[wells], na.rm=T) )
        }
    }
    # Compute the global CV, which is the average of the local ones
    #global_cv = exp(mean(log(replicates_cv[!is.na(replicates_cv)]))) # Geometric mean
    global_cv = mean(replicates_cv, na.rm=TRUE)
    if (is.nan(global_cv) | is.na(global_cv)) { global_cv = default_cv }
    else if (global_cv < min_cv) { global_cv = min_cv }

    # Compute variation for each readout, plus the controls and the blank
    bead_variation = c()
    pvalues = list()
    all_wells = wells_per_treatment
    all_wells[["blank"]] = blanks
    all_wells[["control"]] = controls
    for (treatment in names(all_wells)) {
        wells = all_wells[[treatment]]

        # Wilcox test on replicates
        pvalues[[treatment]] = c()
        if (length(wells) > 1) {
            for ( i in 1:(length(wells)-1) ) {
                for (j in (i+1):length(wells)) {
                    if (length(wilcox_sample[[wells[i]]])!=0 && length(wilcox_sample[[wells[j]]])!=0) {
                        pvalues[[treatment]][paste0(i, "_", j)] = suppressWarnings(wilcox.test(wilcox_sample[[wells[i]]], wilcox_sample[[wells[j]]] + sample(c(-1, 1))/10000, alternative="two.sided")$p.value) # Wilcox test with tiny noise to avoid ties
                    #} else {
                        #pvalues[[treatment]][paste0(i, "_", j)] = 1
                    }
                    #wilcox_antibody[[analyte[which(region==bead)] ]] = c(wilcox_antibody[[ analyte[which(region==bead)] ]], pvalue )
                    #wilcox_cell_line[[lname]] = c(wilcox_cell_line[[lname]], pvalue)
                }
            }
        }

        # Compare the replicate CV with the boostrap CVs and collect the biggest
        square_sum = c()
        for (well in wells) {
            if (!is.nan(median_cv[well]) && !is.na(median_cv[well])) {
                if (global_cv < median_cv[well]) {
                    square_sum = c(square_sum, median_cv[well]^2 * bead_carac[well]^2)
                } else {
                    square_sum = c(square_sum, global_cv^2 * bead_carac[well]^2)
                }
            } else {
                square_sum = 0
            }
        }
        # Take the standard error of the mean and convert it back to a cv (a.k.a cv of the mean)
        for (well in wells) {
            bead_variation[well] = sqrt(sum(square_sum)/length(wells)) / mean(bead_carac[wells])
        }
    }
    return(list(wilcox=pvalues, cv=bead_variation))
}

#' Write lxb data in MIDAS format
#'
#' Write data extracted and processed using readBeads in MIDAS format
#' @param beads Beads data
writeMIDAS <- function(beads, dname, inhibitors, stimulators, midas_dir="midas_data") {
    writeMIDASfile(beads$datas, beads$cv, dname, beads$blanks, beads$controls, beads$wpt, inhibitors, stimulators, midas_dir)
}

#' Write the data in a file in MIDAS format
#'
#' @param datas The values measured. Should a matrix with the wells as rownames and the analyte as column names
#' @param variations The coefficient of variation from the replicates
#' @param dname The name of the data
writeMIDASfile <- function(datas, variations, dname, blanks, controls, wells_per_treatment, inhibitors, stimulators, midas_dir="midas_data") {
    perturbators = c(stimulators, inhibitors)
    if (nchar(midas_dir) == 0) { mitidas_dir="midas_data" }
    if (!grepl("/$", midas_dir)) { midas_dir=paste0(midas_dir, "/") }
    suppressWarnings(dir.create(midas_dir))
    handle = file(paste0(midas_dir, "blunt_", dname, "_MIDAS.csv", sep=""), open="w")
    vhandle = file(paste0(midas_dir, "blunt_", dname, "_MIDAS.var", sep=""), open="w")
    # First line with headers, remove previous file content
    new_line = "ID:type"
    for (treatment in perturbators) {
        new_line = paste(new_line, ",TR:", treatment, sep="")
    }
    new_line = paste(new_line, ",DA:ALL", sep="")
    for (antibody in colnames(datas)) {
        new_line = paste(new_line, ",DV:", antibody, sep="")
    }
    writeLines(new_line, handle)
    writeLines(new_line, vhandle)
    # Blanks
    for (well in blanks) {
        new_line = "blank";
        for (treatment in perturbators) { # TR fields
            new_line = paste(new_line, ",0", sep="")
        }
        new_line = paste(new_line, ",0", sep="") # DA field
        var_line = new_line
        for (antibody in colnames(datas)) { # DV fields
            new_line = paste(new_line, round(datas[well, antibody]), sep=",")
            var_line = paste(var_line, signif(variations[well, antibody], 3), sep=",")
        }
        write(new_line, handle, sep="\n", append=TRUE)
        write(var_line, vhandle, sep="\n", append=TRUE)
    }
    # Controls
    for (well in controls) {
        new_line = "c";
        for (treatment in perturbators) { # TR fields
            new_line = paste(new_line, ",0", sep="")
        }
        new_line = paste(new_line, ",0", sep="") # DA field
        var_line = new_line
        for (antibody in colnames(datas)) { # DV field
            new_line = paste(new_line, round(datas[well, antibody]), sep=",")
            var_line = paste(var_line, signif(variations[well, antibody], 3), sep=",")
        }
        write(new_line, handle, sep="\n", append=TRUE)
        write(var_line, vhandle, sep="\n", append=TRUE)
    }
    # Samples
    for (treatment in names(wells_per_treatment)) {
        if (!(wells_per_treatment[[treatment]] %in% controls || wells_per_treatment[[treatment]] %in% blanks)) {
            # Identify which stimulus and inhibition are used
            stim = stimulators[stimulators %in% unlist(strsplit(treatment, "\\+"))]
            inhib = inhibitors[inhibitors %in% unlist(strsplit(treatment, "\\+"))]
            # Put the perturbation in the MIDAS format
            init_line = "t"
            for (perturbation in perturbators) {
                if (perturbation %in% c(inhib, stim)) {
                    init_line = paste(init_line, "1", sep=",")
                } else {
                    init_line = paste(init_line, "0", sep=",")
                }
            }
            init_line = paste(init_line, "0", sep=",") # DA field

            for (well in wells_per_treatment[[treatment]]) {
                new_line = init_line
                var_line = new_line;
                for (antibody in colnames(datas)) { # DV fields
                    new_line = paste(new_line, round(datas[well, antibody]), sep=",")
                    var_line = paste(var_line, signif(variations[well, antibody], 3), sep=",")
                }         
                write(new_line, handle, sep="\n", append=TRUE)
                write(var_line, vhandle, sep="\n", append=TRUE)
            }
        }
    }

    close(handle)
    close(vhandle)
}

#' Plot the distribution of beads in each well
#'
#' @param region Number of the beads to analyse
#' @param analyte Name of the analytes corresponding to the bead number of region
#' @param tpw List of the treatments applied to each well (names(tpw) is the name of the well)
analyseBeadDistributions <- function(lxb_dataset, region, analyte="", tpw=list()) {
    if (all(analyte == "")) { analyte == region }
    if (length(analyte) != length(region)) { warning("Length of 'analyte' and 'region' differ, setting 'analyte' to 'region'"); analyte=region }
    for (well in names(lxb_dataset)) {
        wdata = lxb_dataset[[well]]
        if (is.null(wdata)) {
            warning(paste0("Well ", well, " does not contain any bead"))
        } else {
            hist_data = lapply( region, function(X) { dd=wdata[wdata[,"Bead ID"]==X,"RP1 Value"]; if (length(dd)==0) {return(0)}; return(dd) })
            names(hist_data) = analyte
            plot(c(-1, 18), rep(0, 2), type="l", col="grey", xlab="log2 Value", main=paste0("Value for well ", well, ifelse(well%in%names(tpw), paste0(" (", tpw[[well]], ")"), "")), xlim=c(0, 17), ylim=c(0, 1.2), ylab="Density")
            tmp=sapply(analyte, function(X) {
                            if (length(hist_data[[X]]) > 20) {
                                lines(density(log2(hist_data[[X]])), col=which(analyte==X))
                            }
                        })
            nbeads = sapply(hist_data, function(X) { ll=length(X); ifelse(is.null(ll), 0, ll)  })
            legend(13, 1.2, paste0(analyte, " (#", nbeads, ")"), col=1:length(analyte), lwd=1)
        }
    }
}

