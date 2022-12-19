# match_rsID takes 4 arguments:
# build: Chromosome build (currently, either "GRCh37" or "GRCh38"). Build is set to "GRCh37", unless specified by the user as "GRCh38"
# input_filename: Name of the file containing the chromosome numbers and variant positions on the chromosomes
# column_chromosome_number: The column number where the data for the chromosome numbers can be found
# column_chromosome_position: The column number where the data for the variant positions can be found

# NB: Make sure BiocManager and biomaRT is installed. If not, uncomment lines 10, 11, and 13, and run.
#_____________________________________________________________________________
#Install BioMart databases if not available
#if (!require("BiocManager", quietly = TRUE))
 #install.packages("BiocManager")

#BiocManager::install("biomaRt")
#_____________________________________________________________________________

match_rsID <- function(build, input_filename, column_chromosome_number, column_chromosome_position){
  #Load bioMart library
  library(biomaRt)
  build = "GRCh37"
  
  if(build == "GRCh37"){
    built_conditions <- c("ENSEMBL_MART_SNP", "hsapiens_snp", "https://grch37.ensembl.org",
                          'refsnp_id', 'chr_name', 'chrom_start', 'allele',
                          'chr_name', 'start', 'end')
  }
  if(build == "GRCh38"){
    built_conditions <- c("ENSEMBL_MART_SNP", "hsapiens_snp", "https://oct2022.archive.ensembl.org",
                          'refsnp_id', 'chr_name', 'chrom_start', 'allele',
                          'chr_name', 'start', 'end')
  }
  
  ## Use the default ENSEMBL datasets
  snp = useMart(biomart = built_conditions[1], 
                   dataset = built_conditions[2], host = built_conditions[3])
  
  print("Loading dataset...", quote = F)
  
  ## Load data into the data frame, SNP_M
  SNP_M <- read.table(input_filename, header = T)
  
  print("Dataset loaded", quote = F)
  print("Fetching rsIDs from ENSEMBL. This might take a while, and will depend on the number of genetic variants...", quote = F)
  
  # Create an empty data frame to save successful searches to
  output_final <- data.frame()
  # Create a data frame for unsuccessful searches
  output_notfound <- data.frame(chromosome = NA, variant_position = NA)
  
  #Query data
  for(i in 1:nrow(SNP_M)){
    
      ## Submit the query
      output_indiv <- getBM(attributes = c(built_conditions[4], built_conditions[5], built_conditions[6], built_conditions[7]),
                      filters = c(built_conditions[8], built_conditions[9], built_conditions[10]), 
                      values = list(SNP_M[i, column_chromosome_number], SNP_M[i, column_chromosome_position], SNP_M[i, column_chromosome_position]), 
                      mart = snp)
    
      # Create a list of chr:positions not found in GRCh37/38 database
      if(nrow(output_indiv) == 0){
        output_notfound[i, 1] <- SNP_M[i, column_chromosome_number]
        output_notfound[i, 2] <- SNP_M[i, column_chromosome_position]
      }
    
      # Add chr:position matches to the data frame, output_final
      else if (nrow(output_indiv) == 1){
        output_final <- rbind(output_final, output_indiv)
      }
    
      # Filter out rsIDs not matching chr:positions
      else if(nrow(output_indiv) > 1){
        j <- 1
        temp <- data.frame(refsnp_id = NA, chr_name = NA, chrom_start = NA, allele = NA)
        while(j <= nrow(output_indiv)){
          if(SNP_M[i,column_chromosome_position] == output_indiv[j, 3]){
            temp <- rbind(temp, list(output_indiv[j, 1], output_indiv[j, 2], output_indiv[j, 3], output_indiv[j, 4]))
            j <- j + 1
          }else{
            j <- j + 1
          }
        }
        
      # Add only chr:position matches to the data frame, output_final
      output_final <- rbind(output_final, temp)
      }
  }
  
  print("Finished fetching rsIDs from ENSEMBL", quote = F)
  print("Removing duplicates...", quote = F)
  
  # Remove duplicates chr:positions from output_final
  output_final <- output_final[!duplicated(output_final), ]
  # Remove any rows containing NA
  output_final <- na.omit(output_final)
  output_notfound <- na.omit(output_notfound)

  print("Writing data to files...", quote = F)
  
  # Create files with the output
  write.table(output_final, "succesful_rsID_matches", row.names = FALSE)
  write.table(output_notfound, "unsuccesful_rsID_matches", row.names = FALSE)
  
  print("Data saved to files", quote = F)
}