
### 13/09/2017

alleleCounts = read.table("/Users/brodriguez/Research/Projects/Pancancer/Germline/Analysis/TEIBA_0.5.7/SourceElements/Frequency/PCAWG/germline_source_elements.alleleCount.tsv", row.names = 1, header = TRUE)
alleleCounts

## Select EUR, ASN and AFR ancestries as we have not enough number of samples for the other groups
targetAncestries <- c("EUR", "ASN", "AFR")
#targetAncestries <- c("EUR", "ASN")

alleleCounts <- alleleCounts[, targetAncestries]
alleleCounts

## Select only rare source elements
#targetSrcElements <- c("1p35.2", "1q23.3", "2q21.3", "3p24.1", "3q26.1", "5q13.1", "7p12.3", "7q31.2", "8p23.1f", "9q22.33", "10q25.1", "11p11.2", "11q14.2", "21q21.1")
#alleleCounts <- alleleCounts[targetSrcElements, ]

# autosomes
EUR_size = 2160 * 2
ASN_size = 437 * 2
AFR_size = 134 * 2
total_chromosomes = c(EUR_size, ASN_size, AFR_size)
#total_chromosomes = c(EUR_size, ASN_size)

smallest_sample_size = min(total_chromosomes)
smallest_sample_size
## Africans are the ones with the minimun sample size: 268

p_values_fisher = c()

# For each source element
for (element in 1:nrow(alleleCounts)){
  
  observed_counts = alleleCounts[element,]
  observed_non_counts = total_chromosomes - observed_counts
  
  print(observed_counts)
  print(observed_non_counts)
  observed_counts_corrected = c()
  
  # For each ancestry
  for (i in 1:ncol(alleleCounts)){
    print(observed_non_counts[i])
    observations = c(rep('L1',observed_counts[i]),rep('no_L1',observed_non_counts[i]))
    counts_sampled = c()
    
    # Make 100 sub-samples with the smallest population size (AFR) computing the allele count each time
    for (j in 1:1000){
      sample = sample(observations,smallest_sample_size)
      counts_sampled = c( counts_sampled , length(which(sample == 'L1')) )
    }
    mean_counts = mean(counts_sampled)
    observed_counts_corrected = c(observed_counts_corrected, mean_counts)
  }
  observed_counts_corrected = round(observed_counts_corrected,0)
  
  ## The expected counts are computed as if the L1 counts would be evenly distributed among the ancestries
  total_observed_counts = sum(observed_counts_corrected)
  number_ancestries = ncol(alleleCounts)
  expected_counts = rep.int(round(total_observed_counts/number_ancestries,0), number_ancestries)

  print("*************")
  print(observed_counts_corrected)
  print(expected_counts)
  print("*************")

  print('observed_counts')
  print(observed_counts)
  print('observed_counts_corrected')
  print(observed_counts_corrected)
  print('expected_counts')
  print(expected_counts)
  
  #Fisher exact test
  print('p_value_fisher')
  contingency_table <- matrix(c(observed_counts_corrected,expected_counts), nrow=2, byrow = T)
  print(contingency_table)

  ## Don??t do the test if all the values in the table are 0
  # set p-value as NA
  if (all(contingency_table == 0)){
    p_values_fisher = c(p_values_fisher,'NA')
  
  ## Do fisher test
  }else{
    f = fisher.test(contingency_table)
    p_values_fisher = c(p_values_fisher,f$p.value)
  }
  print('_______________________________')
}

q_values_fisher = p.adjust(p_values_fisher, method = 'BH')

alleleCounts

final_df = alleleCounts
final_df[ , "p_values_fisher"] <- p_values_fisher
final_df[ , "q_values_fisher"] <- q_values_fisher
final_df

final_df[ , "p_values_fisher"] <- round(as.numeric(as.character(final_df$p_values_fisher)), 3)
final_df[ , "q_values_fisher"] <- round(as.numeric(as.character(final_df$q_values_fisher)), 3)
final_df

write.table( final_df , file='/Users/brodriguez/Research/Projects/Pancancer/Germline/Analysis/TEIBA_0.5.7/SourceElements/Frequency/sourceElements_population_counts_fisher.tsv', quote=FALSE, sep='\t')




