# Using temp files etc created in previous script (code_for_calling_vcf_mtDNA_8SEPT.R)
# Then in R in same directory
setwd("/Users/aleal62p/Dropbox (Otago University)/Michael/vcfs")

# Load library
library(tidyverse)

# Reading in file
temp <- read_tsv("temp",col_names=TRUE)
clades <- read_tsv("sample_mtDNA_clades.txt",col_names=TRUE)

# Ensuring names match up between vcf and clades
names(temp)[c(-1:-9)] <- gsub("-","",gsub("\\.","",gsub("_","",gsub("batch_2_","",names(temp)[c(-1:-9)]))))
clades$ID <- gsub(" ","",(gsub("\\)","",(gsub("\\(","",gsub("-","",gsub("\\.","",gsub("_","",gsub("batch_2_","",clades$ID)))))))))

# This shows that there are some trailing "Rs" for the Verts in the vcf
clades$ID[which(!(clades$ID %in% names(temp)[c(-1:-9)]))]
col_index <- which(grepl("Vert",names(temp)))
names(temp)[col_index] <- gsub("R","",names(temp)[col_index])

# Renaming Takashi
names(temp)[which(names(temp)=="FukurogitsunecleanR")] <- "Takashi"

# Some have a S[0-9][0-9] in front of their code
clades$ID[which(!(clades$ID %in% gsub("S[0-9][0-9]","",names(temp))))]
names(temp) <- gsub("S[0-9][0-9]","",names(temp))

# Also need to filter out a mitogenome sequence from clades that is not RNAseq
clades <- clades %>% filter(ID!="NC003039")

# Need to add entry for Sandy being duplicated
clades %>% filter(ID=="Sandyliver")
clades <- rbind(clades,c("Sandyspleen","RNA_Sandy_spleen", "Blue_NSW"))

# Final double check 
clades$ID[which(!(clades$ID %in% names(temp)[c(-1:-9)]))]
names(temp)[which(!(names(temp) %in% clades$ID))]

# Reference definitions
TAS <- "Takashi"
NSW <- c("Vert35","Vert36","Vert37","Vert38","Vert39","Vert40","Vert41")
KAN <- "SRR3901715"

# Total number of samples
(dim(temp)[2]-9)
# Total number of SNPs
dim(temp)[1]

# Minimum mean depth	Max missing	Total	Mitochondrial genes	Proportion	      Name
# 10    	            5,00%	      38875	7864	              20,2289389067524	Q30DP10mis95_mit.recode.vcf

# For testing out function
#refsample1 <- "Takashi"
#refsample2 <- "Vert35"

# Filter out SNPs that are found on the mitogenome
temp <- temp %>% filter(`#CHROM`!="chrM")

# For testing out function
#refsample1 <- TAS
#refsample2 <- "Vert35"

# For this analysis, only running on TAS vs Vert35 to Vert41 because don't have nuclear reference corresponding to red haplotype

refsample1_refsample2_comparison <- function(refsample1,refsample2) {
  
  output <- matrix(NA, ncol=(dim(temp)[2]-9), nrow=dim(temp)[1])
  
  for (i in 1:dim(temp)[1]) {
    refsample1_genotype <- gsub("\\|","/",gsub(":.*","",temp[i,which(names(temp)==refsample1)]))
    refsample2_genotype <- gsub("\\|","/",gsub(":.*","",temp[i,which(names(temp)==refsample2)]))
    if (refsample1_genotype==refsample2_genotype) {
      # If the genotypes are the same (no resolution)
      output[i,] <- "No_Res"
    } else {
      if (all(c(refsample1_genotype,refsample2_genotype) %in% c("0/0", "1/1"))) {
        # If both references are not hets 
        output[i,(which(gsub("\\|","/",gsub(":.*","",temp[i,]))==refsample1_genotype))-9] <- "ref1"
        output[i,(which(gsub("\\|","/",gsub(":.*","",temp[i,]))==refsample2_genotype))-9] <- "ref2"
        output[i,(which(gsub("\\|","/",gsub(":.*","",temp[i,]))=="0/1"))-9] <- "Het"
      } else { 
        output[i,] <- "No_Res"
      }
    }
  }
  
  output <- gsub("ref1","Green_TAS",output)
  output <- gsub("ref2","Blue_NSW",output)
  
  for (j in 1:dim(output)[1]) {
    if(output[j,1]!="No_Res") {
      for (k in 1:dim(output)[2]) {
        if (output[j,k]!="Het" & !is.na(output[j,k])) {
          if(clades$mtDNA_clade[which(clades$ID==names(temp)[k+9])]==output[j,k]) {
            output[j,k] <- "Match"
          } else {
            output[j,k] <- "Mismatch"
          }
        }
      }
    }
  }
  
  output <- as_tibble(output)
  names(output) <- names(temp)[-c(1:9)]

  output <- as_tibble(cbind(temp[,1:9],output))
  
  output <- output %>% filter(Takashi!="No_Res")
  
  no_samples <- dim(output)[2] - 9
  
  output_total <- output %>% mutate(Match=rowSums(.[10:dim(output)[2]]=="Match",na.rm=TRUE),
                    Mismatch=rowSums(.[10:dim(output)[2]]=="Mismatch",na.rm=TRUE),
                    Het=rowSums(.[10:dim(output)[2]]=="Het",na.rm=TRUE),
                    Total=Match+Mismatch+Het,
                    Prop_geno=Total/no_samples,
                    Match=Match/Total,
                    Mismatch=Mismatch/Total,
                    Het=Het/Total) %>% 
    filter(Het!=0)
  
  remove_cols <- which(names(temp) %in% c("Vert35","Vert36","Vert37","Vert38","Vert39","Vert40","Vert41","Takashi","SRR3901715","SRR8658969","L1","L2","L3","L4"))
  
  output_otago <- output[-remove_cols]
  
  output_otago <- output_otago[-(which(names(output_otago) %in% clades$ID[which(clades$mtDNA_clade=="Red_NA")]))]
  
  no_samples <- dim(output_otago )[2] - 9
  
  output_otago <- output_otago %>% mutate(Match=rowSums(.[10:dim(output_otago)[2]]=="Match",na.rm=TRUE),
                                    Mismatch=rowSums(.[10:dim(output_otago)[2]]=="Mismatch",na.rm=TRUE),
                                    Het=rowSums(.[10:dim(output_otago)[2]]=="Het",na.rm=TRUE),
                                    Total=Match+Mismatch+Het,
                                    Prop_geno=Total/no_samples,
                                    Match=Match/Total,
                                    Mismatch=Mismatch/Total,
                                    Het=Het/Total) %>% 
    filter(Het!=0)
  
  output <- list(output_total,output_otago)
  
  return(output)
  
}    

TAS_Vert35 <- refsample1_refsample2_comparison(TAS,"Vert35")
TAS_Vert36 <- refsample1_refsample2_comparison(TAS,"Vert36")
TAS_Vert37 <- refsample1_refsample2_comparison(TAS,"Vert37")
TAS_Vert38 <- refsample1_refsample2_comparison(TAS,"Vert38")
TAS_Vert39 <- refsample1_refsample2_comparison(TAS,"Vert39")
TAS_Vert40 <- refsample1_refsample2_comparison(TAS,"Vert40")
TAS_Vert41 <- refsample1_refsample2_comparison(TAS,"Vert41")

combined_TAS_Vert <- TAS_Vert35[[2]]
matching_rows <- cbind(TAS_Vert35[[2]][,1:2],"TAS_Vert35[[2]]")

references <- list(TAS_Vert36[[2]],TAS_Vert37[[2]],TAS_Vert38[[2]],TAS_Vert39[[2]],TAS_Vert40[[2]],TAS_Vert41[[2]])
reference_names <- c("TAS_Vert36[[2]]","TAS_Vert37[[2]]","TAS_Vert38[[2]]","TAS_Vert39[[2]]","TAS_Vert40[[2]]","TAS_Vert41[[2]]")

for (k in 1:length(references)) {
  i <- references[[k]]
  for (j in 1:dim(i)[1]) {
    if(nrow(combined_TAS_Vert %>% filter(`#CHROM`==as.character(i[j,1]) & POS==as.numeric(i[j,2])))==0) {
      combined_TAS_Vert <- rbind(combined_TAS_Vert,i[j,])
      matching_rows <- rbind(matching_rows,c(as.matrix(i[j,1:2]),reference_names[k]))
    } else {
      temprow <- combined_TAS_Vert %>% filter(`#CHROM`==as.character(i[j,1]) & POS==as.numeric(i[j,2]))
      if(all(temprow==i[j,],na.rm=TRUE)) {
        matching_rows[(which(matching_rows[,1]==as.character(i[j,1]) & matching_rows[,2]==as.numeric(i[j,2]))),3] <- paste(matching_rows[(which(matching_rows[,1]==as.character(i[j,1]) & matching_rows[,2]==as.numeric(i[j,2]))),3],reference_names[k],sep=",")
      } else {
          print(k,j)
        }
      }
    }
  }

ggplot() + geom_point(data=combined_TAS_Vert,mapping=aes(y=Match,x=Mismatch))
matching_rows <- as_tibble(matching_rows)

# Obvious outlining genes from the graph
combined_TAS_Vert %>% filter(Mismatch<0.2 & Match<0.2) %>%
  select(`#CHROM`,POS,ID,REF,ALT,Match,Mismatch,Het,Prop_geno)
temp %>% filter(`#CHROM`=="chr2" & POS==335113686) %>% as.matrix()
temp %>% filter(`#CHROM`=="chr7" & POS==15766065) %>% as.matrix()
temp %>% filter(`#CHROM`=="chr7" & POS==15767171) %>% as.matrix()

# Genomic location of these
input <- read_tsv("../possum_mito_gene_matches.txt")
input %>% filter(start<335113686 & end>335113686)
input %>% filter(start<15766065 & end>15766065)
input %>% filter(start<15767171 & end>15767171)

# Were they found in more than one individual?
# Yes for this one
matching_rows %>% filter(`#CHROM`=="chr2" & POS==335113686)
# Just with one reference for this one
matching_rows %>% filter(`#CHROM`=="chr7" & POS==15766065)
matching_rows %>% filter(`#CHROM`=="chr7" & POS==15767171)

# Weirdness in variant calling from RNAseq could be explanation
# some of these might need to be verified with genomic sequencing

# Other outlier is one with  a high level of mismatch
combined_TAS_Vert %>% filter(Mismatch>0.65) %>%
  select(`#CHROM`,POS,ID,REF,ALT,Match,Mismatch,Het,Prop_geno)
temp %>% filter(`#CHROM`=="chr1" & POS==39740980) %>% as.matrix()
input %>% filter(start<39740980 & end>39740980)

# We might want to colour each of the dots by the number of refs supporting
# it
names(matching_rows)[3] <- "references"
matching_rows <- matching_rows %>% mutate(POS=as.numeric(POS)) %>% 
                                            arrange(`#CHROM`,POS) %>% 
  mutate(num_refs=nchar(references)-nchar(gsub(",","",references))+1)

combined_TAS_Vert <- combined_TAS_Vert %>%  arrange(`#CHROM`,POS)

combined_TAS_Vert <- full_join(combined_TAS_Vert,matching_rows,by=c("#CHROM","POS"))

# Having a look to see if any patterns are driven by the number of references
ggplot() + geom_point(data=combined_TAS_Vert,mapping=aes(y=Match,x=Mismatch,colour=as.factor(num_refs)))
ggplot() + geom_point(data=(combined_TAS_Vert %>% filter(num_refs>6)),mapping=aes(y=Match,x=Mismatch))

combined_TAS_Vert %>% filter(num_refs>6 & Match>0.5) 
temp %>% filter(`#CHROM`=="chr1" & POS==33713506) %>% as.matrix()
input %>% filter(start<33713506 & end>33713506)

temp %>% filter(`#CHROM`=="chrX" & POS==52282588) %>% as.matrix()
input %>% filter(start<52282588 & end>52282588)

# Making sure the different Sandy samples match
combined_TAS_Vert %>% filter(Sandyliver!=Sandyspleen) %>% select(Sandyliver,Sandyspleen,`#CHROM`,POS)
## FORMAT GT:AD:DP ECT --- DP = APPROXIMATE READ DEPTH ##
## site 1 
temp %>% filter(`#CHROM`=="chr1" & POS==193841171) %>% as.matrix()
#Sandyliver                                   Sandyspleen                                    
# [1,] "0/1:18,7:25:99:.:.:217,0,677:."      "1/1:0,5:5:15:.:.:177,15,0:."
## read depths both below 10 @7 & 5

# site 2
temp %>% filter(`#CHROM`=="chr2" & POS==136442798) %>% as.matrix()
#Sandyliver                          Sandyspleen                      
#  "0/1:14,8:22:99:268,0,519"       "0/0:7,0:7:0:0,0,209"
## spleen has much lower read depth here 

#site 3
temp %>% filter(`#CHROM`=="chr3" & POS==168090302) %>% as.matrix()
# liver                         #spleen
# "0/1:39,25:64:99:852,0,1430"  "0/0:2,0:2:6:0,6,74"
# spleen has read depth of 2 and liver has 99

#site 4
temp %>% filter(`#CHROM`=="chr3" & POS==168090377) %>% as.matrix()
#liver                             spleen
#"0/1:37,18:55:99:586,0,1364"      "0/0:2,0:2:3:0,3,45"
## liver has much higher read depth here 55 - 2 (spleen)

#site 5
temp %>% filter(`#CHROM`=="chr3" & POS==168090452) %>% as.matrix()
#liver                          #spleen
#0/1:52,14:66:99:386,0,1936"    "0/0:2,0:2:6:0,6,84"
# again really low spleen read depths
#site 6
temp %>% filter(`#CHROM`=="chr3" & POS==314659158) %>% as.matrix()
#liver                        #spleen
#"0/1:23,4:27:86:86,0,855"    "0/0:27,0:27:0:0,0,957"
## unclear both have similar read depths here both in favour of reference allele

#site 7
temp %>% filter(`#CHROM`=="chr3" & POS==427745246) %>% as.matrix()
#liver                                      #spleen
# "0/1:58,7,0:65:5:5,0,1465,179,1486,1666"  "0/0:78,8,0:86:48:0,48,2023,235,2047,2233"
#both have read depths above 50 here both have more reads backing the reference allele

#site 8
temp %>% filter(`#CHROM`=="chr4" & POS==257177628) %>% as.matrix()
#liver                    #spleen
#"0/1:4,6:10:99:221,0,136" "0/0:1,0:1:3:0,3,42"
# read depth of 1 here for spleen 10 for liver

#site 9
temp %>% filter(`#CHROM`=="chr4" & POS==257178482) %>% as.matrix()
#liver                       #spleen
#"0/1:12,14:26:99:503,0,419" "1/1:0,2:2:6:84,6,0"
# spleen again another low read depth 

#site 10
temp %>% filter(`#CHROM`=="chr5" & POS==18890489) %>% as.matrix()
#liver                       #spleen
#"0/1:11,9:20:99:312,0,395" "0/0:2,0:2:6:0,6,80"
# another low read depth here for the spleen

#site 11
temp %>% filter(`#CHROM`=="chr5" & POS==41428156) %>% as.matrix()
#liver                       #spleen
#"0/1:23,16:39:99:515,0,836" "1/1:0,2:2:6:84,6,0"
#another low read depth for the spleen

#site 12
temp %>% filter(`#CHROM`=="chr5" & POS==202064268) %>% as.matrix()
#liver                    #spleen
#"0/1:6,3:9:57:57,0,142" "0/0:7,0:9:21:0,21,200"
#both here have a read depth of 9

#site 13
temp %>% filter(`#CHROM`=="chr7" & POS==49270069) %>% as.matrix()
#liver                          #spleen
#"0/1:37,41:78:99:1478,0,1296" "0/0:3,0:3:0:0,0,38"
# low read depth here for spleen

#site 14
temp %>% filter(`#CHROM`=="chr9" & POS==40841564) %>% as.matrix()
#liver                    #spleen
#"0/1:1,3:4:30:113,0,30" "0/0:3,0:3:0:0,0,42"
#low read depth here for both liver and spleen.

## if depth is right after genotype then it seems in at least one of the samples for every site that the read depth is quite low, in some cases 
## both read depths are low 

# Saving table for future use
write.table(combined_TAS_Vert,"combined_match_mismatch.txt",row.names = FALSE,col.names = TRUE,quote=FALSE)

# Saving modified clades for future use
write.table(clades,"sample_mtDNA_clades_renamed.txt",row.names = FALSE,col.names = TRUE,quote=FALSE)
