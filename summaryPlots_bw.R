#############################################
###Function: PLoS_paper_summary.R
###  --------------
###momo.sander@googlemail.com
##############################################

rat_infile <- "../data/interAssay_Rattus_norvegicus_chembl_10.tab"
human_infile <- "../data/interAssay_Homo_sapiens_chembl_10.tab"
para_infile <- "../data/paralogs_chembl_10_compara_62.tab"
ortho_infile <- "../data/orthologs_chembl_10_compara_62.tab"
sample_infile <- "../data/orthologs_chembl_10_compara_62.tab"
compara_paralogs <- "../data/paralogs_62.txt"
compara_orthologs <- "data/orthologs_62.txt"

## @knitr prep
# Loading.
load <- function(infile)
{
   input_table <- read.table(infile ,sep="\t", header=TRUE, na.strings = 'None')
   data <- as.data.frame(input_table) 
   return(data)
}


## @knitr prep_homologs
# Loading.
input_table <- read.table( compara_paralogs ,sep="\t", header=TRUE, na.strings = 'None')
comparaParalogs <- as.data.frame(input_table)

input_table <- read.table( compara_orthologs ,sep="\t", header=TRUE, na.strings = 'None')
comparaOrthologs <- as.data.frame(input_table)

## @knitr filter
# Filtering.
filter <- function(data){
   data$diff <- data$afnty1 - data$afnty2
   data <- data[abs(data$diff) < 20 & data$diff != 0, ]
   return(data)
}

## @knitr dataSummary
# Plot summary data for ortholog and inter-assay comparisons.

density_plot <- function(data, al)
{
   ggplot(data, aes(x=diff))+
   geom_density(fill="darkgrey", col = "darkgrey")+
   scale_x_continuous(breaks=c(-4,-2,0,2,4), limits = c(-5,5))+
   scale_y_continuous(breaks=c(0, 0.2, 0.4,0.6), limits=c(0,.7))+
   geom_rug(col=rgb(.5,0,0,alpha=al))
   
}
scatter_plot <- function(data, al)
{
   ggplot(data, aes(afnty1,afnty2))+
   geom_point(size=1, col = 'darkgrey', shape = 1, alpha = al)+
   scale_x_continuous(breaks=c(2,4,6,8,10,12), limits=c(2,12))+
   scale_y_continuous(breaks=c(2,4,6,8,10,12), limits=c(2,12))+
   coord_equal(ratio = 1)
}


## @knitr bland
bland_plot <- function(data, probs, al)
{
   data$av <- (data$afnty1 + data$afnty2)/2
   ggplot(data)+
   geom_hline(yintercept = quantile(data$diff, prob = probs[3]), col = '#993333')+
   geom_point(aes(x=av, y = diff), col = 'darkgrey', shape =1, size =1, alpha = al)+
   #geom_point(aes(x=av, y = abs(diff)-8), col = 'black', shape = 1, size =1)+
   geom_hline(yintercept = quantile(data$diff, prob = probs[c(2,4)]), lty = 'dotted', col = 'black')+
   geom_hline(yintercept = quantile(data$diff, prob = probs[c(1,5)]), lty = 'dotted', col = 'black')+
   scale_y_continuous(breaks = c(-6,-4,-2,0,2,4,6),limits = c(-6,6))+
   scale_x_continuous(breaks = c(0,2,4,6,8,10,12),limits = c(2,12))
}

## @knitr bootstrap_quantiles
library(boot)
get_q <- function(data,d, probs)
{
   return(quantile(data[d], prob = probs))
}

boot_quant <- function(data, probs, R)
{
   return(boot(data, get_q, probs = probs, R = R))
}


## @knitr take_samples
take_samples <- function(tb, data, intvl, n){
    for (i in intvl) {
        tempFrame <- data[sample(nrow(data), n),]
        tempFrame$it <- i
        tb <- rbind(tb, tempFrame)
    }
    return(tb)
}

## @knitr take_two_samples
take_two_samples <- function(tb, data1, data2, intvl, n){
    for (i in intvl) {
        tempFrame <- rbind(data1[sample(nrow(data1), 1500), ], data2[sample(nrow(data2), 1500), ])
        tempFrame$it <- i
        tb <- rbind(tb, tempFrame)
    }
    return(tb)
}

## @knitr samples
ggplot(tb, aes(x = diff, group = it,  col = it, y = ..density..)) +
geom_density(size = .1, adjust = 1.25) +
coord_cartesian(xlim = c(-5,5))+ 
scale_x_continuous(breaks = c(-4,-2,0,2,4))+
scale_y_continuous(limits = c(0,.65))+
theme(legend.position = "None")

## @knitr summary
ggplot(plotFrame, aes(x = diff,  col = spec, y = ..density..))+
geom_density() +
scale_colour_brewer(palette = "Paired") +
coord_cartesian(ylim = c(0,0.8), xlim = c(-4,4))+
scale_x_continuous(breaks = c(-4,-2,0,2,4))+
theme(legend.position = "None")


ggplot(plotFrame, aes(x = diff,  col = spec, y = ..density..))+ 
geom_density() +
scale_colour_brewer(palette = "Paired") +
coord_cartesian(ylim = c(0,0.8), xlim = c(-4,4))+
scale_x_continuous(breaks = c(-4,-2,0,2,4))+
facet_wrap(~ target_class_L1, ncol=2)+
theme(legend.position = "None")

ggplot(plotFrame, aes(x = target_class_L1,  fill = spec, y = ..count..))+
geom_histogram(position = "dodge")+ 
scale_fill_brewer(palette = "Paired")+
theme(axis.text.x = element_text(angle = 45, hjust = 1))+
theme(legend.position = "None")




## @knitr qq_laplace
library(lawstat)
plot_qq <- function(data)
{
    data$trim <- with(data, cut(diff, breaks=quantile(diff, probs=c(0,1,99,100)/100), labels = c('trim','ok','trim'), include.lowest=TRUE))
    data$ll <- rlaplace(length(data$diff))
    x <- qlaplace(c(.18,.82))
    y <- quantile(data$diff, c(0.18, 0.82))
    slp <- diff(y)/diff(x)
    incp <- mean(data$diff)
    ggplot(data)+
       geom_point(aes(y=sort(diff),x=sort(ll), col = sort(trim)),shape  =1, size =1)+
      geom_abline(slope = slp, lty = 2, col = '#993333', size =.5, intercept = incp)+
       theme(legend.position = 'None')+
       scale_colour_manual(values = c('darkgrey', 'black'))
}


#@knitr qq_normal
plot_qq <- function(data){    
    data$trim <- with(data, cut(diff, breaks=quantile(diff, probs=c(0,1,99,100)/100), labels = c('trim','ok','trim'), include.lowest=TRUE))
    data$nn <- rnorm(data$diff)
    x <- qnorm(c(.2,.8))
    y <- quantile(data$diff, c(0.18, 0.82))
    slp <- diff(y)/diff(x)
    incp <- mean(data$diff)
    ggplot(data)+
       geom_point(aes(y=sort(diff),x=sort(nn), col = sort(trim)), shape = 1, size =1)+
       geom_abline(slope = slp, lty = 2, col = '#993333', size =.5)+
       scale_x_continuous(breaks = c(-8,-6,-4,-2,0,2,4,6,8))+
       scale_y_continuous(breaks = c(8,-6,-4,-2,0,2,4,6,8))+
       theme(legend.position = 'None')+
       scale_colour_manual(values = c('darkgrey', 'black'))
}

    

## @knitr seqDensities
dens_plot <- function(seqDat, diffDat)
{
    data = data.frame(seq_id = seqDat, diff = diffDat)
    ggplot(data,aes(seq_id,abs(diff)))+
    stat_density2d(geom = "tile", aes(fill = ..density..), contour = FALSE, n = 200)+
    coord_cartesian( ylim = c(0,5))+
    scale_fill_gradient(low='white', high = 'black')+
    theme(legend.position = "none")
}

dens_plot(paraFrame$seq_id, paraFrame$diff)
dens_plot(paraFrame$dom_seq_id, paraFrame$diff)
dens_plot(paraFrame$bs_seq_id, paraFrame$diff)
dens_plot(orthoFrame$seq_id, orthoFrame$diff)



## @knitr mwDens
ggplot(paraFrame,aes(molweight,abs(diff)))+
stat_density2d(geom = "tile", aes(fill = ..density..), contour = FALSE, n = 500)+
coord_cartesian( ylim = c(0,5), xlim = c(0,750))+
scale_fill_gradient(low='white', high = 'black')+
theme_bw()



## @knitr seqBox
# Make a SeqId bin label.
binSize <- 20
paraFrame <- paraFrame[order(paraFrame$seq_id),]
paraFrame$seq_id_bin <- sapply(paraFrame$seq_id, function(z) as.integer(z/binSize)*(binSize/10))
paraFrame$dom_seq_id_bin <- sapply(paraFrame$dom_seq_id, function(z) as.integer(z/binSize)*(binSize/10))
paraFrame$bs_seq_id_bin <- sapply(paraFrame$bs_seq_id, function(z) as.integer((z-1)/binSize)*(binSize/10))

ggplot()+
   geom_boxplot(aes(factor(paraFrame$seq_id_bin),abs(paraFrame$diff),
fill = factor(paraFrame$seq_id_bin)),size =1)+
    scale_fill_manual(values=c(rgb(45,0,255,255,maxColorValue=255), rgb(97,0,157,maxColorValue=255),rgb(150,0,104,maxColorValue=255), rgb(202,0,52,maxColorValue=255),rgb(255,0,0,maxColorValue=255)))+
    scale_y_continuous(limits = c(0,5))+
    opts(legend.position = "none")

ggplot()+
   geom_boxplot(aes(factor(paraFrame$dom_seq_id_bin),abs(paraFrame$diff),
fill = factor(paraFrame$dom_seq_id_bin)),size =1)+
    scale_fill_manual(values=c(rgb(45,0,255,255,maxColorValue=255), rgb(97,0,157,maxColorValue=255),rgb(150,0,104,maxColorValue=255), rgb(202,0,52,maxColorValue=255),rgb(255,0,0,maxColorValue=255)))+
    scale_y_continuous(limits = c(0,5))+
    opts(legend.position = "none")

ggplot()+
    geom_boxplot(aes(factor(paraFrame$bs_seq_id_bin),abs(paraFrame$diff),fill = factor(paraFrame$bs_seq_id_bin)),size =1)+
    scale_fill_manual(values=c(rgb(45,0,255,255,maxColorValue=255), rgb(97,0,157,maxColorValue=255),rgb(150,0,104,maxColorValue=255), rgb(202,0,52,maxColorValue=255),rgb(255,0,0,maxColorValue=255),rgb(255,0,0,maxColorValue=255) ))+
    scale_y_continuous(limits = c(0,5))+
    opts(legend.position = "none")



## @knitr mwBox
# Make an Mw bin label.
numBins <- 5
paraFrame <- paraFrame[order(paraFrame$molweight),]
numBreaks <- c(1:(numBins-1))
mwBreaks <- sapply(numBreaks,function(z) paraFrame$molweight[z*length(paraFrame$molweight)/numBins])
paraFrame$mwBin <- findInterval(paraFrame$molweight, mwBreaks)
mwBreaks <- c(0,mwBreaks)
paraFrame$mwBin <- unlist(sapply(paraFrame$mwBin, function(z)  mwBreaks[z+1]))

ggplot()+ 
    geom_boxplot(aes(factor(paraFrame$mwBin),abs(paraFrame$diff)),size =1)+ 
    scale_y_continuous(limits = c(0,5))


## @knitr pairnames
# Aggregate the data for compounds tested against the same pairs.
set_pairnames <- function(data)
{
   data$pairname <- mapply(function(a,b){
          paste(a, b, sep= "_")
        }, as.character(data$accession1), as.character(data$accession2))
   return(data)
}

## @knitr aggregate

agg_frame <- function(data)
{
   data$afnty1Av <- sapply(data$pairname, function(z)
                       median(data$afnty1[data$pairname == z]))
   data$afnty2Av <- sapply(data$pairname, function(z)
                       median(data$afnty2[data$pairname == z]))
   data$numConcat <- sapply( data$pairname, function(z)
         length(data$afnty1[data$pairname == z]))
   aggFrame <- data[!duplicated(data$pairname),]
   return(aggFrame)
}


## @knitr scatterSeq
# Make a scatter plot with dot positions averaged for each pair, col by seqId.
ggplot(aggFrame, aes(afnty1Av,afnty2Av,  size = numConcat, col = seq_id))+
geom_point()+
geom_abline(intercept = 0, slope =1)+
scale_x_continuous(breaks=c(4,6,8,10))+
scale_y_continuous(breaks=c(4,6,8,10))+
scale_colour_gradient(low = 'blue', high = 'red', breaks = c(0,20,40,60,80,100))+
scale_size(trans = 'log')+
coord_cartesian(xlim=c(3.5,10.5), ylim=c(3.5, 10.5))+
coord_equal(ratio = 1)


ggplot(aggFrame, aes(afnty1Av,afnty2Av,  size = numConcat, col = dom_seq_id))+
geom_point()+
geom_abline(intercept = 0, slope =1)+
scale_x_continuous(breaks=c(4,6,8,10))+
scale_y_continuous(breaks=c(4,6,8,10))+
scale_colour_gradient(low = 'blue', high = 'red', breaks = c(0,20,40,60,80,100))+
scale_size(trans = 'log')+
coord_cartesian(xlim=c(3.5,10.5), ylim=c(3.5, 10.5))+
coord_equal(ratio = 1)


ggplot(aggFrame, aes(afnty1Av,afnty2Av,  size = numConcat, col = bs_seq_id))+
geom_point()+
geom_abline(intercept = 0, slope =1)+
scale_x_continuous(breaks=c(4,6,8,10))+
scale_y_continuous(breaks=c(4,6,8,10))+
scale_colour_gradient(low = 'blue', high = 'red', breaks = c(0,20,40,60,80,100))+
scale_size(trans = 'log')+
coord_cartesian(xlim=c(3.5,10.5), ylim=c(3.5, 10.5))+
coord_equal(ratio = 1)


## @knitr scatterTC
# Make a scatter plot with dot positions averaged for each pair, col by TC.
ggplot(aggFrame, aes(afnty1Av,afnty2Av,  size = numConcat, col = target_class_L1))+
   geom_point()+
   geom_abline(intercept = 0, slope =1)+
   scale_x_continuous(breaks=c(2,4,6,8,10,12),  limits=c(4,10))+
   scale_y_continuous(breaks=c(2,4,6,8,10,12),  limits = c(4,10))+
   scale_colour_brewer(palette = 'Dark2')+
   scale_size(trans = 'log')+
   coord_equal(ratio = 1)


## @knitr overlay
# Prepare a data frame for plotting.
plotFrame<- rbind(
data.frame(x="Paralogs in trial, global",seqId=paraFrame$seq_id),
data.frame(x="Orthologs in trial, global",seqId=orthoFrame$seq_id),
data.frame(x="Paralogs in Compara",seqId=comparaParalogs$seqId),           
data.frame(x="Orthologs in Compara",seqId=comparaOrthologs$seqId)
)
ggplot(plotFrame, aes(x=seqId, fill =x))+
   geom_density(aes(x=seqId, y=..density..))+
   scale_fill_manual("",values=c('#003366','#993333', '#003366', '#993333'))+
   theme(legend.position = "none")

sampleVec <- sampleFrame$diff
distFrame <- rbind(
data.frame(x="Paralogs",diff=paraFrame$diff),
data.frame(x="Inter-assay",diff=sampleVec),
data.frame(x="Orthologs",diff=orthoFrame$diff)
)

ggplot(distFrame, aes(x=diff, col = x))+
    geom_density(size = 1)+
    coord_cartesian(xlim = c(-6,6), ylim = c(0,0.6))+
    scale_x_continuous(breaks = c(-6,-3,0,3,6))+
    scale_y_continuous(breaks = c(0,.2,.4, .6))+
    scale_colour_manual(values = c('#003366','darkgrey', '#993333'))+
    theme(legend.position = "none")

## @knitr overlay_paralogs
sampleVec <- sampleFrame$diff
distFrame <- rbind(
data.frame(x="Inter-assay",diff=sampleVec),
data.frame(x="Orthologs",diff=orthoFrame$diff)
)

ggplot(distFrame, aes(x=diff, col = x))+
    geom_density(size = 1)+
    coord_cartesian(xlim = c(-6,6), ylim = c(0,0.6))+
    scale_x_continuous(breaks = c(-6,-3,0,3,6))+
    scale_y_continuous(breaks = c(0,.2,.4, .6))+
    scale_colour_manual(values = c('darkgrey', '#993333'))+
    theme(legend.position = "none")

## @knitr overlay_paralogs
sampleVec <- humanFrame$diff
distFrame <- rbind(
data.frame(x="Paralogs",diff=paraFrame$diff),
data.frame(x="Inter-assay",diff=sampleVec)
)
ggplot(distFrame, aes(x=diff, col = x))+
    geom_density(size = 1)+
    coord_cartesian(xlim = c(-6,6), ylim = c(0,.2,.4,.6))+
    scale_x_continuous(breaks = c(-6,-3,0,3,6))+
    scale_y_continuous(breaks = c(0,.2,.4,.6))+
    scale_colour_manual(values = c('#003366','darkgrey'))+
    theme(legend.position = "none")


## @knitr mle_models
# Maximum likelihood estimation Models.
require(lawstat)
require(fitdistrplus)

mle_laplace <- function(data){
    params <- fitdist(data$diff, "laplace", start=list(scale=1, location = 0))
    laplace<-dlaplace(data$diff, location = as.numeric(params$estimate[2]), scale = as.numeric(params$estimate[1]))
    return(laplace)
}


mle_normal <- function(data){
    params <- fitdistr(data$diff, "normal")
    normal <-dnorm(data$diff, mean = as.numeric(params$estimate[1]), sd = as.numeric(params$estimate[2]))
    return(normal)
}

ortho_laplace <- mle_laplace(orthoFrame)
para_laplace <- mle_laplace(paraFrame)
sample_laplace <- mle_laplace(sampleFrame)

ortho_normal <- mle_normal(orthoFrame)
para_normal <- mle_normal(paraFrame)
sample_normal <- mle_normal(sampleFrame)



## @knitr ana_models
# Analogous Models.
require(lawstat)
require(fitdistr)
ana_laplace <- function(data){
   scale <- mad(data$diff)
   location =  median(orthoFrame$diff)
   laplace<-dlaplace(data$diff, location = location, scale = scale)
   return(laplace)
}

ana_normal <- function(data){
   scale <- sd(data$diff)
   location = mean(orthoFrame$diff)
   normal<-dnorm(data$diff, mean = location, sd = scale)
   return(normal)
}

ortho_laplace <- ana_laplace(orthoFrame)
para_laplace <- ana_laplace(paraFrame)
sample_laplace <- ana_laplace(sampleFrame)

ortho_normal <- ana_normal(orthoFrame)
para_normal <- ana_normal(paraFrame)
sample_normal <- ana_normal(sampleFrame)


## @knitr models
# Plotting models on top of histograms.

model_plot <- function(data, normal, laplace)
{  
   data$normal <- normal
   data$laplace <- laplace
   ggplot(data) + 
   geom_histogram(aes(diff, y= ..density..), fill = "darkgrey", col = "black", binwidth = .4) +
   coord_cartesian(xlim = c(-5, 5), ylim = c(0,.7)) +
   scale_x_continuous(breaks = c(-4,-2,0,2,4))+
   geom_line(aes(x = diff,y = normal), size = .5, col = "#003366")+ #003366
   geom_line(aes(x = diff,y = laplace), size = .5, col = "#993333")
}

model_plot(sampleFrame, sample_normal, sample_laplace)
model_plot(orthoFrame, ortho_normal, ortho_laplace)
model_plot(paraFrame, para_normal, para_laplace)




## @knitr models_dens
# Statistical Models.
library(lawstat)
orthoModel<-dlaplace(orthoFrame$diff, location = 0, scale = 0.85)
paraModel<-dlaplace(paraFrame$diff, location = 0, scale = 1.3)
sampleModel<-dlaplace(sampleFrame$diff, location = 0, scale = 0.85)

model_plot <- function(data, model)
{
    data$model <- model
    ggplot(data)+
    geom_density(aes(diff), fill = "darkgrey", col = "darkgrey")+
    coord_cartesian(xlim = c(-5,5), ylim = c(0,.6))+
    geom_line(aes(x=diff, y =model), size =.5 , col = "#993333")
}
model_plot(sampleFrame, sampleModel)
model_plot(orthoFrame, orthoModel)
model_plot(paraFrame, paraModel)


## @knitr compareMolweight
test_out <- lapply(sort(unique(paraFrame$pairname)), function(z) try(cor.test(paraFrame$molweight[paraFrame$pairname == z], abs(paraFrame$diff[paraFrame$pairname == z]), method = 'spearman'), silent = T))
# Extract and transform p-values.
pValues <- sapply(test_out, function(x) x[3])
pValues <- as.numeric(pValues)
lps <- -log10(pValues)
# Extract and transform yMeans.
yMean <- sapply(test_out, function(x) x[4])
yMean <- as.numeric(yMean)
numSamples <- sapply(sort(unique(paraFrame$pairname)), function(z)
length(sample(paraFrame$diff[paraFrame$pairname == as.character(z)])))
numSamples <- as.numeric(numSamples)
prefName1 <- lapply(sort(unique(paraFrame$pairname)), function(z) as.character(paraFrame$prefName_accession1)[paraFrame$pairname == z][1])
prefName1 <- as.character(prefName1)
prefName2 <- lapply(sort(unique(paraFrame$pairname)), function(z) as.character(paraFrame$prefName_accession2)[paraFrame$pairname == z][1])
prefName2 <- as.character(prefName2)
# Populate the plotFrame.
plotFrame <- data.frame(accession = sort(unique(paraFrame$pairname)),
pValue = pValues, yMean = yMean,  lp = lps, numSamples = numSamples, prefName1 = prefName1, prefName2 = prefName2 )
# Volcano Plot
ggplot()+
   geom_point(aes( x = plotFrame$yMean, y = plotFrame$lp, size = plotFrame$numSamples), col = "darkgrey")+
   scale_size_continuous(name = "# tested")


## @knitr compareTargets
# Calculate p-values.
testResults  <- lapply(sort(unique(wFrame$pairname)), function(z) try(wilcox.test(interSample, wFrame$diff[wFrame$pairname == as.character(z)]
, conf.int = TRUE, exact = FALSE)
        ,silent = TRUE)
)
# Extract and transform p-values.
pValues <- sapply(testResults, function(x) x$p.value)
#pValues <- p.adjust(pValues, method = "bonferroni")
pValues <- as.numeric(pValues)
lps <- -log10(pValues)
# Extract and transform y estimates.
yMean <- sapply(testResults, function(x) x$estimate[[1]])
yMean <- as.numeric(yMean)
yMean <- -yMean
sign <- as.character(lapply(yMean, function(z) isTRUE(z >0)))
# Prepare vectors for plotFrame.
prefName_accession1<- sapply(sort(unique(wFrame$pairname)), function(z) wFrame$prefName_accession1[wFrame$pairname == as.character(z)][1])    
prefName_accession2<- sapply(sort(unique(wFrame$pairname)), function(z) wFrame$prefName_accession2[wFrame$pairname == as.character(z)][1])
seq_id <- sapply(sort(unique(wFrame$pairname)), function(z) wFrame$seq_id[wFrame$pairname == as.character(z)][1])
seq_id <- as.numeric(seq_id)
numSamples <- sapply(sort(unique(wFrame$pairname)), function(z) length(wFrame$diff[wFrame$pairname == as.character(z)]))
numSamples <- as.numeric(numSamples)
# Generate a temporary frame.
plotFrame <- data.frame(pairname = sort(unique(wFrame$pairname)),pValue = pValues, yMean = yMean, sign = sign, numSamples = numSamples, lp = lps, seq_id = seq_id, prefName_accession1 = prefName_accession1, prefName_accession2 = prefName_accession2)

write.table(plotFrame[order(plotFrame$pValue),], file = outfile, sep = '\t', row.names = FALSE, quote=FALSE)
# Volcano Plot
ggplot(plotFrame, aes( x = yMean, y = lp, col = seq_id))+
   geom_point(aes(size = numSamples))+
   #geom_text(aes(label = pairname), col = 'black')+
   coord_cartesian(xlim = c(-4,4))

# Barplot
ggplot(plotFrame, aes(x = sign, y = ..count..))+geom_bar()


## @knitr compareTargetsBootstrap
# Take n samples from the sample distribution.
n = 300
iterations = 1
summFrame <- data.frame(pairname = character(0),pValue = numeric(0), yMean = numeric(0), numSamples = numeric(0), lp = numeric(0), seq_id = numeric(0))
for(i in 1:iterations){
    # Calculate p-values for samples of size n.
    testResults  <- lapply(sort(unique(wFrame$pairname)), function(z) try(wilcox.test(interSample, sample(wFrame$diff[wFrame$pairname == as.character(z)], min(n, length(sample(wFrame$diff[wFrame$pairname == as.character(z)]))))
    , conf.int = TRUE, exact = FALSE)
            ,silent = TRUE)
    )
    # Extract and transform p-values.
    pValues <- sapply(testResults, function(x) x$p.value)
    pValues <- p.adjust(pValues, method = "bonferroni")
    pValues <- as.numeric(pValues)
    lps <- -log10(pValues)
    # Extract and transform y estimates.
    yMean <- sapply(testResults, function(x) x$estimate[[1]])
    yMean <- as.numeric(yMean)
    yMean <- -yMean
    # Prepare vectors for plotFrame.    
    seq_id <- sapply(sort(unique(wFrame$pairname)), function(z) wFrame$seq_id[wFrame$pairname == as.character(z)][1])
    seq_id <- as.numeric(seq_id)
    numSamples <- sapply(sort(unique(wFrame$pairname)), function(z) min(n, length(sample(wFrame$diff[wFrame$pairname == as.character(z)]))))
    numSamples <- as.numeric(numSamples)
    # Generate a temporary frame.
    tmp <- data.frame(pairname = sort(unique(wFrame$pairname)),pValue = pValues, yMean = yMean, numSamples = numSamples, lp = lps, seq_id = seq_id) 
    summFrame <- rbind(summFrame, tmp)
}
write.table(summFrame[order(summFrame$pValue),], file = outfile, sep = '\t', row.names = FALSE, quote=FALSE)

# aggregate into plotFrame
seq_id <- sapply(sort(unique(summFrame$pairname)), function(z)                                       summFrame$seq_id[summFrame$pairname == z][1])
numSamples <- sapply(sort(unique(summFrame$pairname)), function(z)
             summFrame$numSamples[summFrame$pairname == z][1]) 
yMean_md <- sapply(sort(unique(summFrame$pairname)), function(z)                                       median(summFrame$yMean[summFrame$pairname == z]))
yMean_sd <- sapply(sort(unique(summFrame$pairname)), function(z)
                       sd(summFrame$yMean[summFrame$pairname == z]))
lp_md <- sapply(sort(unique(summFrame$pairname)), function(z)
                       median(summFrame$lp[summFrame$pairname == z]))
lp_sd <- sapply(sort(unique(summFrame$pairname)), function(z)
                       sd(summFrame$lp[summFrame$pairname == z]))
plotFrame <- data.frame(pairname = sort(unique(summFrame$pairname)), yMean_md = yMean_md, yMean_sd = yMean_sd, lp_md = lp_md, lp_sd = lp_sd, numSamples = numSamples, seq_id = seq_id)
plotFrame <- subset(plotFrame, !plotFrame$pairname %in% c('P36544_Q05941'))

# Volcano Plot
ggplot(plotFrame, aes( x = yMean_md, y = lp_md, col = seq_id))+
   geom_point(aes(size = numSamples))+
   geom_segment(aes(x=plotFrame$yMean_md-plotFrame$yMean_sd, xend=plotFrame$yMean_md+plotFrame$yMean_sd, yend = lp_md), colour="black")+
   geom_errorbar(aes(ymin=plotFrame$lp_md-plotFrame$lp_sd, ymax=plotFrame$lp_md+plotFrame$lp_sd), colour="black")+
   #geom_text(aes(label = pairname), col = 'black')+
   coord_cartesian(xlim = c(-4,4))



