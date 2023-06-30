# function to calculate Boyce index
# pred = predictions (i.e. 0-1 values from habitat suitability scores)
# res = result = present or absent/background
# bins = number of score bins to use in calculation

BoyceIndex <- function(pred, res, bins=10, bin.method="percentile", range="observed")
{

    if(max(pred)>1 | min(pred)<0)
        {stop("Prediction scores (pred) must be between 0-1")}

    if(sum(res==1) + sum(res==0) != length(res))
        {stop("Presence/Bacground (res) indicators must be 1=presence or 0=background")}

    if(bin.method!="equal" && bin.method!="percentile" )
        {stop("method for binning should be either 'equal' or 'percentile'")}

    if(range!="observed" && range!="expected" )
        {stop("method for determining range should be either 'observed' or 'expected' (0-1)")}

    # set min and max values
    if(range=="expected")
    {
        bi.min<-0
        bi.max<-1
    }
    else #range=="observed"
    {
        bi.min<-min(pred)
        bi.max<-max(pred)
    }

    # find bin thresholds
    if(bin.method=="equal")
    {
        bin.thresholds<-bi.min+(0:bins)/(bins/(bi.max-bi.min))
    }
    else # bin.method="percentile"
    {
        bin.thresholds<-quantile(pred, 0:bins/bins)
    }

    # get bins
    mycuts<-cut(pred,bin.thresholds)

    z<-list("Boyce Index")
    # set up output
    class(z)

    z$method<-bin.method
    z$range<-range

    # set up data frame
    z$bins<-data.frame(row.names=1:bins)
    z$bins[["Bin"]]<-bin.thresholds[1:10]

    # work out total presence and background points
    z$TotalPresences<-sum(res)
    z$TotalBackground<-length(res)-z$TotalPresences

    # count for presences
    bin.pres<-summary(mycuts[res==1])
    bin.abs<-summary(mycuts[res==0])

    # fill counts into output structure
    z$bins[["Observed.N"]]<-bin.pres[1:bins]
    z$bins[["Expected.N"]]<-bin.abs[1:bins]
    # cut function has a bug
    # cut thresholds are rounded, so min value for first bin can exceed
    #  the actual min value, in this case the value is classed n/a
    #  this is then placed in the N+1th bin... so add these to the first bin
    if(length(bin.pres)>bins)
    {
        z$bins[["Observed.N"]][1]<-z$bins[["Observed.N"]][1]+bin.pres[bins+1]
    }
    if(length(bin.abs)>bins)
    {
        z$bins[["Expected.N"]][1]<-z$bins[["Expected.N"]][1]+bin.abs[bins+1]
    }
    z$bins[["Observed.P"]]<-z$bins[["Observed.N"]]/z$TotalPresences
    z$bins[["Expected.P"]]<-z$bins[["Expected.N"]]/z$TotalBackground

    z$bins[["Observed/Expected"]]<-z$bins[["Observed.P"]]/z$bins[["Expected.P"]]

    # finally do correlation of obeserved/expected with rank
    z$BoyceIndex<-cor(z$bins[["Observed/Expected"]], 1:bins, method="spearman")

    z
}
