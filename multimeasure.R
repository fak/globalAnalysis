#############################################
###Function: PLoS_paper_summary.R
###  --------------
###momo.sander@googlemail.com
##############################################

## @knitr sampledFrame
# Loading.
library(ggplot2)
library(MASS)

infile <- '../data/interAssaySampled_human.tab'
intable <- read.table(infile, sep = '\t', header = T)
humanFrame <- as.data.frame(intable)

infile <- '../data/interAssaySampled_rat.tab'
intable <- read.table(infile, sep = '\t', header = T)
ratFrame <- as.data.frame(intable)

#sampledFrame <- rbind(ratFrame[sample(nrow(ratFrame), 1500), ], humanFrame[sample(nrow(humanFrame), 1500), ])
#write.table(sampledFrame,eol = "\n", file = "../data/sampledFrame.tab", sep = '\t', col.names =T, row.names = F, quote = F)
infile <- "../data/sampledFrame.tab"
intable <- read.table(infile, sep = '\t', header = T)
sampledFrame <- as.data.frame(intable)


## @knitr interFrame
inter_infile <- "../data/inter_assay_pairs_00.tab"
input_table <- read.table(inter_infile, sep = "\t", header = TRUE, na.strings = "None")
interFrame <- as.data.frame(input_table)


## @knitr filter
sampledFrame$diff <- sampledFrame$afnty1 - sampledFrame$afnty2
sampledFrame <- sampledFrame[abs(sampledFrame$diff) < 20 & sampledFrame$diff != 0, ]

interFrame$diff <- interFrame$afnty1 - interFrame$afnty2
interFrame <- interFrame[abs(interFrame$diff) < 20 & interFrame$diff != 0, ]

## @knitr compare
sampled <- sampledFrame$diff
median <- interFrame$diff
median_ex <- interFrame$diff[interFrame$measured == 2]
distFrame <- rbind(data.frame(x = "median", diff = interFrame$diff), data.frame(x = "median excluded", diff = median_ex), data.frame(x = "sampled", diff = sampled))

ggplot(distFrame, aes(x = diff, col = x, y = ..density..)) + geom_density(size = .5) + coord_cartesian(xlim = c(-5,5)) + scale_colour_manual(values = c("lightblue", "darkblue", "darkgrey"))+theme_bw()

ggplot(distFrame, aes(x = diff, fill = x, y = ..density..)) + geom_histogram( position = 'dodge') + coord_cartesian(xlim = c(-5,5)) + scale_fill_manual(values = c("lightblue", "darkblue", "darkgrey"))+theme_bw()



## @knitr bootstrap_sampled
iterations <- 100
sampledFrame$it <- iterations * 10
tb <- sampledFrame
for(i in 1:iterations){

   tempFrame <- rbind(ratFrame[sample(nrow(ratFrame), 1500), ],humanFrame[sample(nrow(humanFrame), 1500), ])

   tempFrame$diff <- tempFrame$afnty1 - tempFrame$afnty2
   tempFrame <- tempFrame[abs(tempFrame$diff) < 20 & tempFrame$diff != 0, ]
   tempFrame$it <- i
   tb <- rbind(tb, tempFrame)
}

ggplot(tb, aes(x = diff, group = it,  col = it, y = ..density..)) + geom_density(size = .1) + coord_cartesian(xlim = c(-5,5))+theme_bw()



## @knitr bootstrap_inter
iterations <- 100
interFrame$it <- iterations * 10
tb <- interFrame
tb$tid <- NULL
tb$target_class_L1 <- "dummy"
for(i in 1:iterations){

   tempFrame <- rbind(ratFrame[sample(nrow(ratFrame), 1500), ],humanFrame[sample(nrow(humanFrame), 1500), ])

   tempFrame$diff <- tempFrame$afnty1 - tempFrame$afnty2
   tempFrame <- tempFrame[abs(tempFrame$diff) < 20 & tempFrame$diff != 0, ]
   tempFrame$it <- i
   tempFrame$uniprot <- NULL
   tb <- rbind(tb, tempFrame)
}

ggplot(tb, aes(x = diff, group = it,  col = it, y = ..density..)) + geom_density(size = .1) + coord_cartesian(xlim = c(-5,5))+theme_bw()





 

