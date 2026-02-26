compute_bsa <- function(pdb_path = pdb_path, domain_map = NULL, cleanup = TRUE){
  
  pdb <- read.pdb(pdb_path)
  chainA <- atom.select(pdb, chain = "A")
  chainB <- atom.select(pdb, chain = "B")
  pdbA <- trim.pdb(pdb, chainA)
  pdbB <- trim.pdb(pdb, chainB)
  
  write.pdb(pdbA, file = "A.pdb")
  write.pdb(pdbB, file = "B.pdb")
  write.pdb(pdb,  file = "complex.pdb")
  
  # Run FreeSASA with RSA output
  system("freesasa A.pdb --format=rsa --probe-radius=1.4 -o A.rsa")
  system("freesasa B.pdb --format=rsa --probe-radius=1.4 -o B.rsa")
  system("freesasa complex.pdb --format=rsa --probe-radius=1.4 -o complex.rsa")
  
  # Parse RSA files
  parse_rsa <- function(file) {
    lines <- readLines(file)
    res_lines <- grep("^RES", lines, value = TRUE)
    fields <- do.call(rbind, strsplit(res_lines, "\\s+"))
    data.frame(
      chain = fields[, 2],
      resn  = fields[, 3],
      resno = as.integer(fields[, 4]),
      sasa  = as.numeric(fields[, 5]),
      stringsAsFactors = FALSE
    )
  }
  
  
  rsaA <- parse_rsa("A.rsa")
  rsaB <- parse_rsa("B.rsa")
  rsaC <- parse_rsa("complex.rsa")
  
  # Combine and subtract for per-residue BSA
  combine_bsa <- function(monomer, complex) {
    merged <- merge(monomer, complex, by = c("chain", "resn", "resno"), suffixes = c(".mono", ".complex"))
    merged$bsa <- merged$sasa.mono - merged$sasa.complex
    merged[merged$bsa > 0, c("chain", "resn", "resno", "bsa")]
  }
  
  bsaA <- combine_bsa(rsaA, rsaC)
  bsaB <- combine_bsa(rsaB, rsaC)
  per_residue_bsa <- rbind(bsaA, bsaB)
  
  tmp <- rbind(rsaA, rsaB)
  tmp$sasa <- 0
  for(ii in 1:nrow(per_residue_bsa)){
    ind <- intersect(x = which(tmp$resn == per_residue_bsa$resn[ii]), 
                     y = which(tmp$resno == per_residue_bsa$resno[ii]))
    tmp$sasa[ind] <- per_residue_bsa$bsa[ii]
  }
  colnames(tmp) <- c("chain", "resn", "resno", "bsa")
  per_residue_bsa <- tmp
  
  if((is.null(domain_map)) || (nrow(per_residue_bsa) == 0)){
    
    return(per_residue_bsa)
    
  } else {
    
    domains <- rep("PF00000", nrow(per_residue_bsa))
    domain_map <- domain_map[complete.cases(domain_map), ]
    if(nrow(domain_map) > 0){
      for(ii in 1:nrow(per_residue_bsa)){
        
        for(jj in 1:nrow(domain_map)){
          
          if(per_residue_bsa$resn[ii] == domain_map$chain[jj]){
            
            if(per_residue_bsa$resno[ii] %in% (domain_map$pfam_start[jj]:domain_map$pfam_end[jj])){
              
              domains[ii] <- domain_map$pfam[jj]
              
            }
            
          }
          
        }
        
      }
    }
    
    per_residue_bsa$domains <- domains
    
    return(per_residue_bsa)
    
  }
  
  if(cleanup){
    
    file.remove(c("A.pdb", "A.rsa", "B.pdb", "B.rsa",  "complex.pdb", "complex.rsa"))
    
  }
  
}