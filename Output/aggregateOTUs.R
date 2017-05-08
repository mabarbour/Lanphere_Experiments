#how to aggregate OTUs by taxonomy. for example, if you want to see how bacteria respond
#to treatment at the Phylum or Order level.
#make sure metagenomeSeq is loaded
library('metagenomeSeq')

#run through previous script to get a filtered, normalized otu table, 'hb'

hbPhylum <- aggTax(hb, lvl = Phylum, out = "matrix")

#can also aggregate by sample
aggSam()
