#!/usr/bin/env Rscript
library(ggplot2)
library(plyr)
library("RColorBrewer", lib.loc="/vol/r/3.2")

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandArgs(trailingOnly=FALSE))

inputfile <- args[1]

print(inputfile)

funnydata <- read.table(inputfile,  sep="\t", header=TRUE)


###funnydata <- read.table("/home/tamal/PycharmProjects/crappycrap1/overlappingmarkers_summary_egg",  sep="\t", header=TRUE)

col_headings <- c('quorum','indels','found')

col_names_mid <- c('quorum','indels','ref_clust_id',"hit_50.")
col_names_mida <- c('quorum','indels','ref_clust_id',"hit_70.")
col_names_midb <- c('quorum','indels','ref_clust_id',"hit_100.")

funnydata1 <- aggregate(funnydata$"hit_50." ~ funnydata$quorum + funnydata$indels + funnydata$ref_clust_id , data=funnydata, FUN = max )
names(funnydata1) <- col_names_mid
funnydata2 <- aggregate(funnydata1$"hit_50." ~ funnydata1$quorum + funnydata1$indels  , data=funnydata1, FUN = sum )
names(funnydata2) <- col_headings



bp <- ggplot(funnydata2, aes(funnydata2$indels, funnydata2$quorum , fill = funnydata2$found, label = sprintf("%d", funnydata2$found))) + 
  geom_raster(hjust=0.5,vjust=0.5, interpolate=FALSE) + geom_text() + scale_fill_gradientn(colors = brewer.pal(30, "YlGnBu"),na.value = "gray",limits=c(0,29)) +
  xlab("indels") + ylab("quorum") + ggtitle("Cluster findings 50%") + scale_x_continuous(breaks=seq(0,9,1)) + scale_y_continuous(breaks=seq(0,13,1))

bp

ggsave("Clusterfindings50_style2.pdf", width=6, height=6)

dev.off()


 funnydata1a <- aggregate(funnydata$"hit_70." ~ funnydata$quorum + funnydata$indels + funnydata$ref_clust_id , data=funnydata, FUN = max )
 names(funnydata1a) <- col_names_mida
 funnydata2a <- aggregate(funnydata1a$"hit_70." ~ funnydata1a$quorum + funnydata1a$indels  , data=funnydata1a, FUN = sum )
 names(funnydata2a) <- col_headings

 bp <- ggplot(funnydata2a, aes(funnydata2a$indels, funnydata2a$quorum , fill = funnydata2a$found, label = sprintf("%d", funnydata2a$found))) +
   geom_raster(hjust=0.5,vjust=0.5, interpolate=FALSE) + geom_text() + scale_fill_gradientn(colors = brewer.pal(30, "YlGnBu"),na.value = "gray",limits=c(0,29)) +
   xlab("indels") + ylab("quorum") + ggtitle("Cluster findings 70%") + scale_x_continuous(breaks=seq(0,9,1)) + scale_y_continuous(breaks=seq(0,13,1))

 bp
 ggsave("Clusterfindings70_style2.pdf", width=6, height=6)

 dev.off()


 funnydata1b <- aggregate(funnydata$"hit_100." ~ funnydata$quorum + funnydata$indels + funnydata$ref_clust_id , data=funnydata, FUN = max )
 names(funnydata1b) <- col_names_midb
 funnydata2b <- aggregate(funnydata1b$"hit_100." ~ funnydata1b$quorum + funnydata1b$indels  , data=funnydata1b, FUN = sum )
 names(funnydata2b) <- col_headings

 bp <- ggplot(funnydata2b, aes(funnydata2b$indels, funnydata2b$quorum , fill = funnydata2b$found, label = sprintf("%d", funnydata2b$found))) +
   geom_raster(hjust=0.5,vjust=0.5, interpolate=FALSE) + geom_text() + scale_fill_gradientn(colors = brewer.pal(30, "YlGnBu"),na.value = "gray",limits=c(0,29)) +
   xlab("indels") + ylab("quorum") + ggtitle("Cluster findings 100%") + scale_x_continuous(breaks=seq(0,9,1)) + scale_y_continuous(breaks=seq(0,13,1))

 bp

 ggsave("Clusterfindings100_style2.pdf", width=6, height=6)

 dev.off()

