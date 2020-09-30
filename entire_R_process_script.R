
###read in sample names file so can use for import commands later
samples_df <- read.csv("samples.csv")


###major library for commands involved in vcf analysis
library(VariantAnnotation)


###combining sample names and data folder filepath then importing vcfs
tmp_filename <- paste0("data/", samples_df[1,]$filename)

tmp_vcf <- readVcf(tmp_filename)


###exploration of data appearance and format
###header(tmp_vcf)
###samples(header(tmp_vcf))
###geno(header(tmp_vcf))
###head(rowRanges(tmp_vcf), 3)
###names(head(rowRanges(tmp_vcf), 3))

###testing how to split up the names in vcf used to refer to loci, doesn't use chrI etc
unique(t(as.data.frame(strsplit(names(rowRanges(tmp_vcf)), ":")))[,1])

###making a data frame called seqnames that will hold the format vcf uses and the name we want
###names such as chrI and chrII
chr <- paste0("chr", as.roman(1:16))
chr <- c(chr, "chrmt")

seqnames_df <- data.frame(ref = unique(seqnames(tmp_vcf)),
                          chr = chr)

###write that dataframe as a csv
write.csv(seqnames_df, "seqnames.csv")

###create a list and then put vcf files into it using readvcf of each iterated filename
vcf_list <- list()
for (id in samples_df$sample) {
  filename <- paste0("data/", subset(samples_df, sample == id)$filename)
  vcf_list[[id]] <- readVcf(filename)
}

###check if variant numbers make sense for strain background (they do for YJM789)
vcf_list
lapply(vcf_list, length)

###filter list for PASS and check for the lengths again
vcf_list_filtered <- lapply(vcf_list, subset, FILTER == "PASS")
lapply(vcf_list_filtered, length)

###create a dataframe showing the effect of filtering on the variant numbers for the vcfs
tmp1_df <- t(as.data.frame(lapply(vcf_list, length)))
tmp2_df <- t(as.data.frame(lapply(vcf_list_filtered, length)))

filter_effect_df <- cbind(tmp1_df, tmp2_df)
colnames(filter_effect_df) <- c("before", "after")

###write this filter effect data frame to a csv for reference later
write.csv(filter_effect_df, "filter_effect.csv")

###Create a list for variants by strain. Original analysis plan was for multiple parent strains
###This is less important for analysis with a single parent strain like this one.
###However this step is before finding the shared variants that are due to being YJM789
###Important because those variants are just due to strain mapping to S288c reference, not in lab evolution.
variants_id_per_strain_list <- list()
for (strain_name in levels(samples_df$strain)) {
  variants_id_per_strain_list[[strain_name]] <- list()
  for (sample in subset(samples_df, strain == strain_name)$sample) {
    variants_id_per_strain_list[[strain_name]][[sample]] <- names(vcf_list_filtered[[sample]])
  }
}

###print a summary of list by strain
for (strain_name in names(variants_id_per_strain_list)) {
  cat(strain_name, ":\n")
  print(summary(variants_id_per_strain_list[[strain_name]]))
  cat("\n")
}

###identify the intersection of variants that we can remove based on strain, not experiment
strain_specific_variants_list <- list()
for (strain_name in names(variants_id_per_strain_list)) {
  strain_specific_variants_list[[strain_name]] <- Reduce(intersect, variants_id_per_strain_list[[strain_name]])
}

###looks like we will be removing 57154 variants...shared between all 16? seems like not enough
###can we remove variants that are found in the lab SNP list?
summary(strain_specific_variants_list)

###create a list of variants that removes all the strain specific
relevant_vcf_list <- list()
for (strain_name in names(strain_specific_variants_list)) {
  for (id in subset(samples_df, strain == strain_name)$sample) {
    if (id %in% names(vcf_list_filtered)) {
      variants_to_keep <- setdiff(names(vcf_list_filtered[[id]]),
                                  strain_specific_variants_list[[strain_name]])
      relevant_vcf_list[[id]] <- vcf_list_filtered[[id]][variants_to_keep]
    }
  }
}

###create a data frame with variant counts from all lists so far
tmp3_df <- data.frame(total = t(as.data.frame(lapply(vcf_list, length))),
                      sample = rownames(t(as.data.frame(lapply(vcf_list, length))))
)

tmp4_df <- data.frame(filtered = t(as.data.frame(lapply(vcf_list_filtered, length))),
                      sample = rownames(t(as.data.frame(lapply(vcf_list_filtered, length))))
)

tmp5_df <- data.frame(relevant = t(as.data.frame(lapply(relevant_vcf_list, length))),
                      sample = rownames(t(as.data.frame(lapply(relevant_vcf_list, length))))
)

variant_count_df <- merge(tmp3_df, tmp4_df, by="sample")
variant_count_df <- merge(variant_count_df, tmp5_df, by="sample")

###write that data frame to a csv
write.csv(variant_count_df, "variant_count.csv")

###making a data frame of relevant variants from relevant vcf
###H1 refers to a sample omitted from another analysis
###There is no H1 in this analysis list, so logic still applies as desired.
###In effect, no reason to alter the logic of the function.
relevant_variants_df_list <- list()
for (id in samples_df$sample) {
  if (id != "H1") {
    variant_ids <- sort(names(relevant_vcf_list[[id]]))
    relevant_variants_df_list[[id]] <- data.frame(tmp1 = variant_ids,
                                                  tmp2 = rep(1, length(variant_ids)))
    colnames(relevant_variants_df_list[[id]]) <- c("variant",id)        
  }    
}

###see what this list looks like
summary(relevant_variants_df_list)

for (id in names(relevant_variants_df_list)) {
  cat(id, ":\n")
  print(summary(relevant_variants_df_list[[id]]))
  print(nrow(relevant_variants_df_list[[id]]))
  cat("\n")
}


###make a matrix of all the variants and their presence or lack in each sample
relevant_variants_df <- Reduce(function(x,y) merge(x = x, y = y, by = "variant", all=TRUE), 
                               relevant_variants_df_list)

relevant_variants_df[is.na(relevant_variants_df)] <- 0

relevant_variants_df

###data fram manipulation library
library(dplyr)

rownames(relevant_variants_df) <- relevant_variants_df$variant
relevant_variants_df <- relevant_variants_df %>% select(-variant)
relevant_variants_df

###improve sorting of variants
relevant_variants <- rownames(relevant_variants_df)
tmp6_df <- t(as.data.frame(strsplit(relevant_variants, ":")))

seq_name <- tmp6_df[,1]
head(seq_name)
head(tmp6_df)
tmp7_df <- t(as.data.frame(strsplit(tmp6_df[,2], "_")))
head(tmp7_df)
tmp7_df[,1] <- as.numeric(tmp7_df[,1])
head(tmp7_df)
coor_start <- as.numeric(tmp7_df[,1])
head(coor_start)


###make a data frame that can be sorted by the coor_start info
variants_info_df <- data.frame(variant = relevant_variants,
                               seq_name = seq_name,
                               start = coor_start)
###check progress
head(variants_info_df)
tail(variants_info_df)

###sort the dataframe
variants_info_df <- variants_info_df %>% arrange(seq_name, start)

###check sort
head(variants_info_df)
tail(variants_info_df)

###start using seqnames_df to rename rows for chromosome number

rownames(seqnames_df) <- seqnames_df$ref
chr <- seqnames_df[variants_info_df$seq_name, ]$chr
###check that all levels of chr numbers have been added
head(chr)
tail(chr)
###add matched chr number column to sorting dataframe which has refnumbers
variants_info_df$chr <- chr
head(variants_info_df)
tail(variants_info_df)

###write this sorting dataframe to a csv to save
write.csv(variants_info_df, "variants_info.csv")

###make the relevant variants matrix sorted by sorting dataframe
relevant_variants_df <- relevant_variants_df[as.character(variants_info_df$variant),]
#relevant_variants_df
#tail(head(relevant_variants_df))

###write the sorted relevant variants dataframe to a csv
write.csv(relevant_variants_df, "relevant_variants.csv")


###clear out some environmental variables
rm(chr, coor_start, filename, id, ref, sample, seq_name, tmp_df, tmp_filename,  tmp_list, tmp_var, tmp_vcf, tmp1_df, tmp2_df, tmp3_df)
rm(var1, var2, var3, variant_ids, variants_to_keep)
rm(tmp4_df, tmp5_df, tmp6_df, tmp7_df)



###PCA ANALYSES###
library(pcadapt)

###make matrix from the relevant variants
relevant_variants_matrix <- read.pcadapt(relevant_variants_df, type = "pcadapt")



relevant_variants_res <- pcadapt(relevant_variants_matrix, K = 15, ploidy = 1)

###previous command gave error, said could not compute SVD
###asked if any SNPs or individuals had only missing values
###used following code to check if any rows = 0
###none did so I tried lower K value for the previous command
###is it a memory issue? because K=15 worked and values beyond 15 do not work.
for(x in 1:7112) {
  if (sum(relevant_variants_df[x,]) == 0) {
    print(x)
  }
}

###screeplot for previous to check for K
plot(relevant_variants_res, option = "screeplot")

###trying K=3
relevant_variants_res <- pcadapt(relevant_variants_matrix, K = 3)
###K=3 gave appropriate screeplot


library(dplyr)

###check for usage of sample_df in plot series names
samples_df$sample == colnames(relevant_variants_df)

###PC1 and 2 plot labeled with strains as resistant/control
plot(relevant_variants_res, 
     option = "scores", 
     pop = samples_df$resistance)
###same plot labeled with strains as YPD/MCHM
plot(relevant_variants_res, 
     option = "scores", 
     pop = samples_df$medium)


###PC1 and PC3
plot(relevant_variants_res, 
     option = "scores",
     i = 1,
     j = 3,
     pop = samples_df$resistance)

###PC2 and PC3
plot(relevant_variants_res, 
     option = "scores",
     i = 2,
     j = 3,
     pop = samples_df$resistance)

###changed pop = samples_df$filename from resistance to make a few more plots to check weirdness


###do LD check even though manhattan plot isn't obviously weird
par(mfrow = c(1, 2))
for (i in 1:2)
  plot(relevant_variants_res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

###same K=20 error as before, used K=15 to create this LD clumping corrected matrix
relevant_variants_no_ld_res <- pcadapt(relevant_variants_matrix, K = 15, LD.clumping = list(size = 200, thr = 0.1), ploidy = 1)
plot(relevant_variants_no_ld_res, option = "screeplot")


###used K=3 again to make the res, it looks best
relevant_variants_no_ld_res <- pcadapt(relevant_variants_matrix, K = 3, LD.clumping = list(size = 200, thr = 0.1), ploidy = 1)


###manhattan plot of this new no ld version
plot(relevant_variants_no_ld_res , option = "manhattan")

###score plots of new no ld version
par(mfrow = c(1, 2))
for (i in 1:2)
  plot(relevant_variants_no_ld_res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

###new plots by resistance for no ld version
plot(relevant_variants_no_ld_res, 
     option = "scores", 
     pop = samples_df$resistance)

plot(relevant_variants_no_ld_res, 
     option = "scores",
     i = 1,
     j = 3,
     pop = samples_df$resistance)

plot(relevant_variants_no_ld_res, 
     option = "scores",
     i = 2,
     j = 3,
     pop = samples_df$resistance)


###Analysis by media

###MCHM
relevant_variants_MCHM_df <- relevant_variants_df %>% select(as.character(subset(samples_df, medium=="MCHM")$sample))


relevant_variants_MCHM_matrix <- read.pcadapt(relevant_variants_MCHM_df, type = "pcadapt")

###K=7 was highest that returned no error
relevant_variants_MCHM_res <- pcadapt(relevant_variants_MCHM_matrix, K = 7, ploidy = 1)

###K=2? Need to ask Amaury
plot(relevant_variants_MCHM_res, option = "screeplot")

relevant_variants_MCHM_res <- pcadapt(relevant_variants_MCHM_matrix, K = 2, ploidy = 1)

plot(relevant_variants_MCHM_res, 
     option = "scores", 
     pop = subset(samples_df, medium=="MCHM")$filename)

###YPD
relevant_variants_YPD_df <- relevant_variants_df %>% select(as.character(subset(samples_df, medium=="YPD")$sample))


relevant_variants_YPD_matrix <- read.pcadapt(relevant_variants_YPD_df, type = "pcadapt")

###K=7 was highest that returned no error
relevant_variants_YPD_res <- pcadapt(relevant_variants_YPD_matrix, K = 7, ploidy = 1)

###K=2
plot(relevant_variants_YPD_res, option = "screeplot")

relevant_variants_YPD_res <- pcadapt(relevant_variants_YPD_matrix, K = 2, ploidy = 1)

plot(relevant_variants_YPD_res, 
     option = "scores", 
     pop = subset(samples_df, medium=="YPD")$filename)



###code for keeping rows from the dataframe that only have rowsums above or below a number
###allows pca analysis limited to variants that are shared by a number of strains,
###excluding those that are unique to only one strain, etc. First analysis looks at
###variants that are in 1 at minimum and 5 strains at most. Goal to analyze
###pca to check if patterns of clustering by resistance pop up if we remove variants shared
###by too many strains, ie masking patterns of variants that might be causing resistance in
###resistant strains only. Nothing popped up. The orignal PCA of all 6800 variants was most informative.
relevant_variants_under5_df <- relevant_variants_df[rowSums(relevant_variants_df) <= 5,]


relevant_variants_under5_matrix <- read.pcadapt(relevant_variants_under5_df, type = "pcadapt")

###still same K=15 maximum
relevant_variants_under5_res <- pcadapt(relevant_variants_under5_matrix, K = 15, ploidy = 1)

###K=2
plot(relevant_variants_under5_res, option = "screeplot")


relevant_variants_under5_res <- pcadapt(relevant_variants_under5_matrix, K = 2, ploidy = 1)

plot(relevant_variants_under5_res, 
     option = "scores", 
     pop = samples_df$resistance)


###code for keeping rows from the dataframe that only have rowsums above or below a number
relevant_variants_under3_df <- relevant_variants_df[rowSums(relevant_variants_df) <= 3,]


relevant_variants_under3_matrix <- read.pcadapt(relevant_variants_under3_df, type = "pcadapt")

###still same K=15 maximum
relevant_variants_under3_res <- pcadapt(relevant_variants_under3_matrix, K = 15, ploidy = 1)

###K=3
plot(relevant_variants_under3_res, option = "screeplot")


relevant_variants_under3_res <- pcadapt(relevant_variants_under3_matrix, K = 3, ploidy = 1)

plot(relevant_variants_under3_res, 
     option = "scores", 
     pop = samples_df$resistance)

plot(relevant_variants_under3_res, 
     option = "scores",
     i = 1,
     j = 3,
     pop = samples_df$resistance)

plot(relevant_variants_under3_res, 
     option = "scores",
     i = 2,
     j = 3,
     pop = samples_df$resistance)



###code for keeping rows from the dataframe that only have rowsums above or below a number
relevant_variants_under10_df <- relevant_variants_df[rowSums(relevant_variants_df) <= 10,]


relevant_variants_under10_matrix <- read.pcadapt(relevant_variants_under10_df, type = "pcadapt")

###still same K=15 maximum
relevant_variants_under10_res <- pcadapt(relevant_variants_under10_matrix, K = 15, ploidy = 1)

###K=2?
plot(relevant_variants_under10_res, option = "screeplot")


relevant_variants_under10_res <- pcadapt(relevant_variants_under10_matrix, K = 2, ploidy = 1)


plot(relevant_variants_under10_res, 
     option = "scores", 
     pop = samples_df$resistance)


###try a version with only variants shared by at least 2 so no unique variants and under 10 so no strainwide
relevant_variants_2to10_df <- relevant_variants_under10_df[rowSums(relevant_variants_under10_df) >= 2,]


relevant_variants_2to10_matrix <- read.pcadapt(relevant_variants_2to10_df, type = "pcadapt")


relevant_variants_2to10_res <- pcadapt(relevant_variants_2to10_matrix, K = 15, ploidy = 1)


plot(relevant_variants_2to10_res, option = "screeplot")

###K=3?
relevant_variants_2to10_res <- pcadapt(relevant_variants_2to10_matrix, K = 3, ploidy = 1)


plot(relevant_variants_2to10_res, 
     option = "scores", 
     pop = samples_df$resistance)



relevant_variants_5to10_df <- relevant_variants_under10_df[rowSums(relevant_variants_under10_df) >= 5,]


relevant_variants_5to10_matrix <- read.pcadapt(relevant_variants_5to10_df, type = "pcadapt")


relevant_variants_5to10_res <- pcadapt(relevant_variants_5to10_matrix, K = 15, ploidy = 1)


plot(relevant_variants_5to10_res, option = "screeplot")

###K=2?
relevant_variants_5to10_res <- pcadapt(relevant_variants_5to10_matrix, K = 2, ploidy = 1)


plot(relevant_variants_5to10_res, 
     option = "scores", 
     pop = samples_df$resistance)



###variants shared by exactly 8 strains
relevant_variants_exactly8_df <- relevant_variants_df[rowSums(relevant_variants_df) == 8,]


relevant_variants_exactly8_matrix <- read.pcadapt(relevant_variants_exactly8_df, type = "pcadapt")


relevant_variants_exactly8_res <- pcadapt(relevant_variants_exactly8_matrix, K = 15, ploidy = 1)


plot(relevant_variants_exactly8_res, option = "screeplot")

###K=3?
relevant_variants_exactly8_res <- pcadapt(relevant_variants_exactly8_matrix, K = 3, ploidy = 1)


plot(relevant_variants_exactly8_res, 
     option = "scores", 
     pop = samples_df$resistance)


plot(relevant_variants_exactly8_res, 
     option = "scores",
     i = 1,
     j = 3,
     pop = samples_df$resistance)

plot(relevant_variants_exactly8_res, 
     option = "scores",
     i = 2,
     j = 3,
     pop = samples_df$resistance)




####Manipulate relevant_variants_df to remove all variants that are in YJM789 SNPs list
####This is a previous list of known SNPs between YJM789 and the S288c reference genome
####This should catch any remaining SNP variants that did not get filtered out earlier.

YJM789_SNPs <- read.csv("YJM789_SNPs.csv")


###recreating some dataframes for sorting by chr name instead of ref|etc|
#tmp_df <- t(as.data.frame(strsplit(relevant_variants, ":")))
#tmp1_df <- t(as.data.frame(strsplit(tmp_df[,2], "_")))
#tmp1_df[,1] <- as.numeric(tmp1_df[,1])
#coor_start <- as.numeric(tmp1_df[,1])
#variants_info_df <- data.frame(variant = relevant_variants,
#                               seq_name = seq_name,
#                               start = coor_start)

###made relevant variants data frame with chr number and start coordinate for intersecting with YJM789 SNPs
relevant_variants_df[,17] <- variants_info_df$chr
relevant_variants_df[,18] <- variants_info_df$start
relevant_variants_names_df <- relevant_variants_df[, c(17,18,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]

###combine names of chromosomes and start coordinates into a single column so can find intersection by that column
relevant_variants_names_df[,1] <- paste(relevant_variants_names_df[,1], relevant_variants_names_df[,2])
YJM789_SNPs[,1] <- paste(YJM789_SNPs[,1], YJM789_SNPs[,2])
###identify SNPs want to remove because are just YJM789 SNPs
relevant_variants_YJMSNPs <- merge(relevant_variants_names_df, YJM789_SNPs, by.x = 1, by.y = 1, all = FALSE)

###final data frame with only the variants not in the merged dataframe above
final_relevant_variants <- relevant_variants_names_df[!(relevant_variants_names_df$V17 %in% relevant_variants_YJMSNPs$V17),]

write.csv(final_relevant_variants, "final_relevant_variants.csv")

###cannot edit final_relevant_variants easily with variants_inf_df anymore because rows no longer line up
###some were removed


###quick PCA plot of the new final_relevant_variants data frame, must remove first two columns to use as matrix
final_relevant_variants_forpca <- final_relevant_variants[,-1]
final_relevant_variants_forpca <- final_relevant_variants_forpca[,-1]
final_relevant_matrix <- read.pcadapt(final_relevant_variants_forpca, type = "pcadapt")


final_relevant_res <- pcadapt(final_relevant_matrix, K = 15, ploidy = 1)


plot(final_relevant_res, option = "screeplot")

###K=3
final_relevant_res <- pcadapt(final_relevant_matrix, K = 3, ploidy = 1)


plot(final_relevant_res, 
     option = "scores", 
     pop = samples_df$resistance)

plot(final_relevant_res, 
     option = "scores",
     i = 1,
     j = 3,
     pop = samples_df$resistance)

plot(final_relevant_res, 
     option = "scores",
     i = 2,
     j = 3,
     pop = samples_df$resistance)

###same results which is not surprising, only 238 variants were removed using YJM789 SNP list###
###This is the pca analysis that is included in the figures for the paper


###code template to filter out any variants from a dataframe based on one of the strains having it or not
#RM11_outliers_resistance_df <- RM11_outliers_resistance_df %>% 
#  rownames_to_column('variant') %>% 
#  filter(F9==0)  %>%
#  column_to_rownames('variant')

####Code for the rest of the page was executed to produce
####csv files per MCHM ILE strain that only had variants not in controls
####and did coding analysis on them.
####This code was also executed to find variants in the control strains
####and do the coding analysis on them as well to produce csvs for each
####control strain. The only part excluded was the next portion of filtering S1-S8 out
####of final_relevant_variants dataframe. Instead, coding analysis of control
####strains picks up at extracting variants out of final_relevant_variants for each control strain using
###S1_variants <- final_relevant_variants %>% 
###rownames_to_column('variant') %>% 
###filter(S1==1)  %>%
###column_to_rownames('variant')
####then following to
###vcf_S1_variants <- relevant_vcf_list[["S1"]][rownames(S1_variants)]
####The remaining code executes with proper name replacement with S1-S8 variable names 


library(tibble)

###make a new dataframe to remove control snps, making one with only variants in the 8 MCHM ILEs
variants_not_in_controls <- final_relevant_variants %>% 
  rownames_to_column('variant') %>% 
  filter(S1==0)  %>%
  column_to_rownames('variant')

###now iterate through by removing variants in S2-S8
variants_not_in_controls <- variants_not_in_controls %>% 
  rownames_to_column('variant') %>% 
  filter(S2==0)  %>%
  column_to_rownames('variant')

###did all that all the way up to S8
###write variants_not_in_controls to csv
write.csv(variants_not_in_controls, "variants_not_in_controls.csv")


###add a column to add up how many evolved strains had that variant using rowsum
###remove name columns again because they can't be rowsummed
variants_not_in_controls_rowsum <- variants_not_in_controls[,-1]
variants_not_in_controls_rowsum <- variants_not_in_controls_rowsum[,-1]
variants_not_in_controls_rowsum[,17] <- rowSums(variants_not_in_controls_rowsum)
colnames(variants_not_in_controls_rowsum)[17] <- "Sum"
write.csv(variants_not_in_controls_rowsum, "variants_not_in_controls_rowsum.csv")

###remove the chr name columns from variants_not_in_controls to prepare for extracting vcfs from vcf list
###this call will be based on rownames on a dataframe list with only relevant variants for that strain on it
###making this dataframe then making 8 individual dataframes that remove all the rows not in that strain
final_variants_to_make_strain_vcfs <- variants_not_in_controls[,-1]
final_variants_to_make_strain_vcfs <- final_variants_to_make_strain_vcfs[,-1]


S9_variants_not_in_controls <- final_variants_to_make_strain_vcfs %>% 
  rownames_to_column('variant') %>% 
  filter(S9==1)  %>%
  column_to_rownames('variant')

S10_variants_not_in_controls <- final_variants_to_make_strain_vcfs %>% 
  rownames_to_column('variant') %>% 
  filter(S10==1)  %>%
  column_to_rownames('variant')

S11_variants_not_in_controls <- final_variants_to_make_strain_vcfs %>% 
  rownames_to_column('variant') %>% 
  filter(S11==1)  %>%
  column_to_rownames('variant')

S12_variants_not_in_controls <- final_variants_to_make_strain_vcfs %>% 
  rownames_to_column('variant') %>% 
  filter(S12==1)  %>%
  column_to_rownames('variant')

S13_variants_not_in_controls <- final_variants_to_make_strain_vcfs %>% 
  rownames_to_column('variant') %>% 
  filter(S13==1)  %>%
  column_to_rownames('variant')

S14_variants_not_in_controls <- final_variants_to_make_strain_vcfs %>% 
  rownames_to_column('variant') %>% 
  filter(S14==1)  %>%
  column_to_rownames('variant')

S15_variants_not_in_controls <- final_variants_to_make_strain_vcfs %>% 
  rownames_to_column('variant') %>% 
  filter(S15==1)  %>%
  column_to_rownames('variant')

S16_variants_not_in_controls <- final_variants_to_make_strain_vcfs %>% 
  rownames_to_column('variant') %>% 
  filter(S16==1)  %>%
  column_to_rownames('variant')

write.csv(S9_variants_not_in_controls, "S9_variants_not_in_controls.csv")
write.csv(S10_variants_not_in_controls, "S10_variants_not_in_controls.csv")
write.csv(S11_variants_not_in_controls, "S11_variants_not_in_controls.csv")
write.csv(S12_variants_not_in_controls, "S12_variants_not_in_controls.csv")
write.csv(S13_variants_not_in_controls, "S13_variants_not_in_controls.csv")
write.csv(S14_variants_not_in_controls, "S14_variants_not_in_controls.csv")
write.csv(S15_variants_not_in_controls, "S15_variants_not_in_controls.csv")
write.csv(S16_variants_not_in_controls, "S16_variants_not_in_controls.csv")

###extract vcf for each MCHM strain using above dataframes to reference rows to extract
vcf_S9_variants_not_in_controls <- relevant_vcf_list[["S9"]][rownames(S9_variants_not_in_controls)]

vcf_S10_variants_not_in_controls <- relevant_vcf_list[["S10"]][rownames(S10_variants_not_in_controls)]

vcf_S11_variants_not_in_controls <- relevant_vcf_list[["S11"]][rownames(S11_variants_not_in_controls)]

vcf_S12_variants_not_in_controls <- relevant_vcf_list[["S12"]][rownames(S12_variants_not_in_controls)]

vcf_S13_variants_not_in_controls <- relevant_vcf_list[["S13"]][rownames(S13_variants_not_in_controls)]

vcf_S14_variants_not_in_controls <- relevant_vcf_list[["S14"]][rownames(S14_variants_not_in_controls)]

vcf_S15_variants_not_in_controls <- relevant_vcf_list[["S15"]][rownames(S15_variants_not_in_controls)]

vcf_S16_variants_not_in_controls <- relevant_vcf_list[["S16"]][rownames(S16_variants_not_in_controls)]


###go through and change the seqlevel names to chromosome names from ref numbers
seqlevels(vcf_S9_variants_not_in_controls) <- as.character(seqnames_df$chr)
seqlevels(vcf_S10_variants_not_in_controls) <- as.character(seqnames_df$chr)
seqlevels(vcf_S11_variants_not_in_controls) <- as.character(seqnames_df$chr)
seqlevels(vcf_S12_variants_not_in_controls) <- as.character(seqnames_df$chr)
seqlevels(vcf_S13_variants_not_in_controls) <- as.character(seqnames_df$chr)
seqlevels(vcf_S14_variants_not_in_controls) <- as.character(seqnames_df$chr)
seqlevels(vcf_S15_variants_not_in_controls) <- as.character(seqnames_df$chr)
seqlevels(vcf_S16_variants_not_in_controls) <- as.character(seqnames_df$chr)

writeVcf(vcf_S9_variants_not_in_controls, "S9_variants_not_in_controls.vcf")
writeVcf(vcf_S10_variants_not_in_controls, "S10_variants_not_in_controls.vcf")
writeVcf(vcf_S11_variants_not_in_controls, "S11_variants_not_in_controls.vcf")
writeVcf(vcf_S12_variants_not_in_controls, "S12_variants_not_in_controls.vcf")
writeVcf(vcf_S13_variants_not_in_controls, "S13_variants_not_in_controls.vcf")
writeVcf(vcf_S14_variants_not_in_controls, "S14_variants_not_in_controls.vcf")
writeVcf(vcf_S15_variants_not_in_controls, "S15_variants_not_in_controls.vcf")
writeVcf(vcf_S16_variants_not_in_controls, "S16_variants_not_in_controls.vcf")


####Coding Analysis of Variants####
###load transcript library object.
###This version is sacCer3 genome (2011) at UCSC not R64-2-1 (2015)
###6692 vs 6691 transcript rows, 7034 vs 7052 exon rows
###Also loaded the UCSC sacCer3 genome instead of R62-2-1
###release date is 2011 instead of 2014

library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

library(BSgenome.Scerevisiae.UCSC.sacCer3)

s288c_genome <- BSgenome.Scerevisiae.UCSC.sacCer3

###make genomes of vcfs match my genome objects
genome(vcf_S9_variants_not_in_controls) <- genome(txdb)
genome(vcf_S10_variants_not_in_controls) <- genome(txdb)
genome(vcf_S11_variants_not_in_controls) <- genome(txdb)
genome(vcf_S12_variants_not_in_controls) <- genome(txdb)
genome(vcf_S13_variants_not_in_controls) <- genome(txdb)
genome(vcf_S14_variants_not_in_controls) <- genome(txdb)
genome(vcf_S15_variants_not_in_controls) <- genome(txdb)
genome(vcf_S16_variants_not_in_controls) <- genome(txdb)

###predict coding effects
###throwing error related to subsets of seqnames...so trying to get everything to match those in txdb which cannot be easily edited
###in particular seqnames in txdb use chrM instead of chrmt
###edit seqnames2_df to be seqnames_df but with chrM instead of chrmt
seqnames2_df <- seqnames_df
seqnames2_df$chr <- factor(seqnames2_df$chr, levels = c(levels(seqnames2_df$chr), "chrM"))
seqnames2_df[17,2] <- "chrM"

seqnames(s288c_genome) <- as.character(seqnames2_df$chr)
seqlevels(vcf_S9_variants_not_in_controls) <- as.character(seqnames2_df$chr)
seqlevels(vcf_S10_variants_not_in_controls) <- as.character(seqnames2_df$chr)
seqlevels(vcf_S11_variants_not_in_controls) <- as.character(seqnames2_df$chr)
seqlevels(vcf_S12_variants_not_in_controls) <- as.character(seqnames2_df$chr)
seqlevels(vcf_S13_variants_not_in_controls) <- as.character(seqnames2_df$chr)
seqlevels(vcf_S14_variants_not_in_controls) <- as.character(seqnames2_df$chr)
seqlevels(vcf_S15_variants_not_in_controls) <- as.character(seqnames2_df$chr)
seqlevels(vcf_S16_variants_not_in_controls) <- as.character(seqnames2_df$chr)

###actually predict the coding effects now that the error is fixed
###So apparently there are warnings about some of the variants being out of ranges
###Worried this means the genomes releases are different enough that you can't use them
###the loci being used don't match????? Therefore all the annotations are wrong too????
coding_vcf_S9_variants_not_in_controls <- predictCoding(vcf_S9_variants_not_in_controls, txdb, seqSource=s288c_genome)
coding_vcf_S10_variants_not_in_controls <- predictCoding(vcf_S10_variants_not_in_controls, txdb, seqSource=s288c_genome)
coding_vcf_S11_variants_not_in_controls <- predictCoding(vcf_S11_variants_not_in_controls, txdb, seqSource=s288c_genome)
coding_vcf_S12_variants_not_in_controls <- predictCoding(vcf_S12_variants_not_in_controls, txdb, seqSource=s288c_genome)
coding_vcf_S13_variants_not_in_controls <- predictCoding(vcf_S13_variants_not_in_controls, txdb, seqSource=s288c_genome)
coding_vcf_S14_variants_not_in_controls <- predictCoding(vcf_S14_variants_not_in_controls, txdb, seqSource=s288c_genome)
coding_vcf_S15_variants_not_in_controls <- predictCoding(vcf_S15_variants_not_in_controls, txdb, seqSource=s288c_genome)
coding_vcf_S16_variants_not_in_controls <- predictCoding(vcf_S16_variants_not_in_controls, txdb, seqSource=s288c_genome)

###fix this column to be as.character?
coding_vcf_S9_variants_not_in_controls$ALT <- lapply(coding_vcf_S9_variants_not_in_controls$ALT, as.character)
coding_vcf_S10_variants_not_in_controls$ALT <- lapply(coding_vcf_S10_variants_not_in_controls$ALT, as.character)
coding_vcf_S11_variants_not_in_controls$ALT <- lapply(coding_vcf_S11_variants_not_in_controls$ALT, as.character)
coding_vcf_S12_variants_not_in_controls$ALT <- lapply(coding_vcf_S12_variants_not_in_controls$ALT, as.character)
coding_vcf_S13_variants_not_in_controls$ALT <- lapply(coding_vcf_S13_variants_not_in_controls$ALT, as.character)
coding_vcf_S14_variants_not_in_controls$ALT <- lapply(coding_vcf_S14_variants_not_in_controls$ALT, as.character)
coding_vcf_S15_variants_not_in_controls$ALT <- lapply(coding_vcf_S15_variants_not_in_controls$ALT, as.character)
coding_vcf_S16_variants_not_in_controls$ALT <- lapply(coding_vcf_S16_variants_not_in_controls$ALT, as.character)

###write them to tables
write.table(coding_vcf_S9_variants_not_in_controls, "coding_vcf_S9_variants_not_in_controls.csv", sep=";", row.names = FALSE)
write.table(coding_vcf_S10_variants_not_in_controls, "coding_vcf_S10_variants_not_in_controls.csv", sep=";", row.names = FALSE)
write.table(coding_vcf_S11_variants_not_in_controls, "coding_vcf_S11_variants_not_in_controls.csv", sep=";", row.names = FALSE)
write.table(coding_vcf_S12_variants_not_in_controls, "coding_vcf_S12_variants_not_in_controls.csv", sep=";", row.names = FALSE)
write.table(coding_vcf_S13_variants_not_in_controls, "coding_vcf_S13_variants_not_in_controls.csv", sep=";", row.names = FALSE)
write.table(coding_vcf_S14_variants_not_in_controls, "coding_vcf_S14_variants_not_in_controls.csv", sep=";", row.names = FALSE)
write.table(coding_vcf_S15_variants_not_in_controls, "coding_vcf_S15_variants_not_in_controls.csv", sep=";", row.names = FALSE)
write.table(coding_vcf_S16_variants_not_in_controls, "coding_vcf_S16_variants_not_in_controls.csv", sep=";", row.names = FALSE)



