#code for making normalized OTU tables
#load metagenomeSeq package
library(metagenomeSeq)

#samples = columns = a willow plant (ie, 9UF3)
#otus (aka features) = rows = reads that are a 97% match to each other. 
#the matrices (otu tables) exported with this code will only have OTU names like
# "OTU_1", but can be exported to include taxonomy. 

#load fungi otu table with taxonomy
hf<- load_biom("Data/otu_table_wTax_blast.biom")

#load bacteria otu table with taxonomy
hb <- load_biom("Data/otu_table_wTax16s.biom")

#before normalization, filter samples and otus by depth of coverage
kingdom <- (fData(hf)$taxonomy1 == "k__Fungi") #remove OTUS that didn't have a BLAST match to fungi
hf <- hf[kingdom, ]

hf <- filterData(hf, present = 1, depth = 9000)# filters  OTUS present in at least
#one sample (considering removing these singletons) and keeps samples with a read count of at least 9000. 
#I determined these numbers via data exploration and it means we lose 3 samples. 

#saving the filtered OTU table
hffiltmat <- MRcounts(hf, norm=FALSE, log=FALSE)
write.csv(hffiltmat, "hfFilteredCounts.csv")


#filtering OTUs that matched bacteria & archaea, 
#keeping samples with 6000 reads or more (removes one sample).
unique(fData(hb)$taxonomy1)
kingdom <- (fData(hb)$taxonomy1 == "Unclassified")
kingdom <- kingdom == "FALSE"
hb <- hb[kingdom, ] #there are only a few

hb <- filterData(hb, present = 1, depth = 6000)#bactera, similar rationale 
### IMPORTANT!!! before doing any further analysis, I would remove bacteria OTU_1
# and OTU_2. These OTUs are chloroplasts and mitochondria, respectively, and 
# it doesn't make sense to include these in downstream analyses.
# OTU_1 = row 1, OTU_2 = row 4
hb<-hb[-c(1,4), ]

#saving filtered OTU table
hbfiltmat <- MRcounts(hb, norm=FALSE, log=FALSE)
write.csv(hbfiltmat, "hbFilteredCounts.csv")

#this function normalizes OTU tables. 
#depth of sampling varies by orders of magnitude, people who actually care about
#these statistics have developed a statistically appropriate method, called
#cumulative sum scaling, to reduce the effect of sequencing depth /uneven sampling
#on downstream analyses
hf <- cumNorm(hf)
hb <- cumNorm(hb)

#fungi normalized matrix, log-transformed
hfMat<- MRcounts(hf, norm=TRUE, log=TRUE)
hist(rowSums(hfMat), main = "Histogram of normalized fungal OTU abundance")
hist(colSums(hfMat), main = "Histogram of normalized fungal sample size")
head(hfMat)
range(colSums(hfMat))
write.csv(hfMat, "hfNormalizedCounts.csv")

#bateria normalized matrix, log-transformed
hbMat<- MRcounts(hb, norm=TRUE, log=TRUE)
hist(rowSums(hbMat), main = "Histogram of normalized bacterial OTU abundance")
hist(colSums(hbMat), main = "Histogram of normalized bacterial sample size")
write.csv(hbMat, "hbNormalizedCounts.csv")



