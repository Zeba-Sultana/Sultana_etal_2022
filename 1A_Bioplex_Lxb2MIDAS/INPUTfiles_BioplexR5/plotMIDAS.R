#!/usr/local/bin/Rscript
#-*-encoding: utf8 -*-

library("STASNet")

if (!exists("cargs")) {
    cargs = commandArgs(TRUE)
}
fixed_range = FALSE
lim = Inf
show_values = TRUE
for (argument in cargs) {
    if (grepl("^-l", argument)) {
        fixed_range = TRUE
        lim = abs(as.numeric(gsub("^-l", "", argument)))
    } else if (grepl("--nv", argument)) {
        show_values = FALSE
    } else if (grepl(".csv$", argument)) {
        file_name = paste0(argument)
        data_name = gsub(".csv$", "", basename(argument))
        data_name = gsub("_MIDAS", "", data_name)

        full_datas=STASNet::readMIDAS(file_name)
        control=full_datas[grep("^c$|control", rownames(full_datas)),]
        if (!is.null(nrow(control)) && nrow(control) > 1) {
            control = colMeans(control)
        }
        blank=full_datas[rownames(full_datas)=="blank",]
        datas=full_datas[-c( which(rownames(full_datas)=="blank"), grep("control|^c$", rownames(full_datas)) ),]
        control = matrix(rep(control, nrow(datas)), nrow=nrow(datas), byrow=T)
        rep_variation = aggregate(full_datas, by=list(rownames(full_datas)), sd)
        rownames(rep_variation) = rep_variation[,1]
        rep_variation = rep_variation[,-1]
        rep_mean = aggregate(full_datas, by=list(rownames(full_datas)), mean)[,-1]

        pdf(paste0("heatmap_", data_name, ".pdf"), width=5+ncol(datas)/3, height=4+nrow(datas)/6)
        STASNet:::plotHeatmap(log2(datas/control), main=paste0("Log-fold change ", data_name), lim=lim, fixedRange=fixed_range, show_values=show_values)
        full_datas[full_datas==0] = NaN
        STASNet:::plotHeatmap(log2(full_datas), main=paste0("Raw data ", data_name), fixedRange=TRUE, lim=10)
        if (sum(is.na(rep_variation)) < 0.5 * nrow(rep_variation) * ncol(rep_variation)) {
            STASNet:::plotHeatmap(rep_variation/rep_mean, main=paste0("Replicate variation ", data_name))
        }
        dev.off()
    } else if (grepl(".var$", argument)) {
        file_name = paste0(argument)
        data_name = gsub(".var$", "", basename(argument))
        data_name = gsub("_MIDAS", "", data_name)

        
        datas=STASNet::readMIDAS(file_name)
        pdf(paste0("error_heatmap_", data_name, ".pdf"), width=5+ncol(datas)/3, height=4+nrow(datas)/6)
        STASNet:::plotHeatmap(datas, main=paste0("Error ", data_name), col=colorRampPalette(c("white", "red1")))
        dev.off()
    }
}
