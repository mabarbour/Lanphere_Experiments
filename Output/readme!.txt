file explanations, sorry there are so many

hfNormalizedCounts.csv --> use this to get started. hf = humboldt fungi. counts have been filtered and normalized
	as described in the R code: Code4normalizedOTUtables

hbNormalizedCounts.csv --> normalized bacteria otu table. (same process as above)

otus tables(matrices) are "species" matrices...sequence reads are grouped at 97% similarity to define otus/'species'

otu_table_wTax16S.biom --> a .biom file that can be opened in R with various packages, 
	R code in folder shows how this file was used to make normalized otu tables. 
	this file is raw otu data and contains taxonomy data for the otus. 
	16S = reads are from bacteria 16S subunit, aka bacteria otus


otu_table_wTax_blast.biom --> same as above, but _blast indicates fungal ITS sequence data, aka fungal otus

otu_table_ITS --> raw fungal otu data in csv format. same data as biom file but no taxonomy.

otu_table_16S --> ibid, bacterial data


