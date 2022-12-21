match_rsID <- function(build, input_filename, column_chromosome_number, column_chromosome_position){
  #Load bioMart library
  library(biomaRt)
  
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
  
  #Query data
  for(i in 1:nrow(SNP_M)){
    
      ## Submit the query
      output <- getBM(attributes = c(built_conditions[4], built_conditions[5], built_conditions[6], built_conditions[7]),
                filters = c(built_conditions[8], built_conditions[9], built_conditions[10]), 
                values = list(SNP_M[i, column_chromosome_number], SNP_M[i, column_chromosome_position], SNP_M[i, column_chromosome_position]), 
                mart = snp)
    
      # Create a list of chr:positions not found in GRCh37/38 database, and immediately write to file using append
      if(nrow(output) == 0){
        output_notfound <- data.frame(chromosome = NA, variant_position = NA)
        output_notfound[i, 1] <- SNP_M[i, column_chromosome_number]
        output_notfound[i, 2] <- SNP_M[i, column_chromosome_position]
        output_notfound <- na.omit(output_notfound)
        write.table(output_notfound, "unsuccesful_rsID_matches.txt", append = T, row.names = F, col.names = F)
      }
    
      # Immediately write rsID chr:position matches to file using append
      else if (nrow(output) == 1){
        write.table(output, "succesful_rsID_matches.txt", append = T, row.names = F, col.names = F)
      }
    
      # Filter out rsIDs not matching chr:positions and immediately write rsID chr:position matches to file using append
      else if(nrow(output) > 1){
        j <- 1
        temp <- data.frame(refsnp_id = NA, chr_name = NA, chrom_start = NA, allele = NA)
        while(j <= nrow(output)){
          if(SNP_M[i,column_chromosome_position] == output[j, 3]){
            temp <- rbind(temp, list(output[j, 1], output[j, 2], output[j, 3], output[j, 4]))
            j <- j + 1
          }else{
            j <- j + 1
          }
        }
        
      # Add only rsID chr:position matches to the file
      write.table(output, "succesful_rsID_matches.txt", append = T, row.names = F, col.names = F)
      }
  }
  
  print("Finished fetching rsIDs from ENSEMBL", quote = F)
  print("Adding headers and removing duplicates...", quote = F)

  success <- read.table("succesful_rsID_matches.txt", header = F)
  colnames(success) <- c("refsnp_id", "chr_name", "chrom_start", "alleles")
  success<- success[!duplicated(success), ]
  success <- na.omit(success)
  write.table(success, "succesful_rsID_matches.txt", row.names = F)
  fail <- read.table("unsuccesful_rsID_matches.txt", header = F)
  colnames(fail) <- c("chr_name", "chrom_start")
  fail <- fail[!duplicated(fail), ]
  fail <- na.omit(fail)
  write.table(fail, "unsuccesful_rsID_matches.txt", row.names = F)
  
  print("Done!", quote = F)
}