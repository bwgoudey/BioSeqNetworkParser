classify_acc <- function(acc) {
  
  #Regexs are derived from a range of sources
  # https://community.gep.wustl.edu/repository/course_materials_WU/annotation/Genbank_Accessions.pdf
  # https://www.uniprot.org/help/accession_numbers
  # PIR are cooked up. 
  # https://registry.identifiers.org/registry/insdc
  # https://registry.identifiers.org/registry/pdb
  id_regexs=list(
    uniprot_protein="^(([A-N,R-Z][0-9]([A-Z][A-Z, 0-9][A-Z, 0-9][0-9]){1,2})|([O,P,Q][0-9][A-Z, 0-9][A-Z, 0-9][A-Z, 0-9][0-9]))(.\\d+)?$",
    genbank_protein="^[A-Z]{3}(\\d{5}|\\d{7})(\\.\\d+)?$",
    genbank_nucleotide="^(([UJMA]\\d{5})|([A-Z]{2}\\d{6}))(\\.\\d+)?$",
    genbank_contig="^(([A-Z]{4}\\d{8}(\\d+)?)|([A-Z]{6}\\d{9}))(\\.\\d+)?$",
    refseq_protein="^((AP|NP|XP|YP|ZP)_\\d+)(\\.\\d+)?$",
    refseq.anon_protein="^((WP)_\\d+)(.\\d+)?$",
    refseq_nucleotide="^(((AC|NC|NG|NM|NR|NT|NW|XM|XR)_\\d+)|(NZ\\_[A-Z]{2,4}\\d+))(\\.\\d+)?$",
    pdb_protein="^[0-9][A-Za-z0-9]{3}(_([ABCLX]?))?$"#,
    #pir_protein="^(([A-Z]{2}\\d{4})|([A-OR-TV-Z]{1}\\d{5}))$"
  )
  
  #Setup a  dataframe containing the return values
  # Base it directly on names of regex bank to ensure no mapping issues
  vals=data.frame(r=names(id_regexs)) %>% separate(r, sep="_", into=c("db", "seq_type"))
  df=data.frame(db=rep(NA,length(acc)), seq_type=rep(NA,length(acc)))
  #Match the regex
  I = which(sapply(id_regexs, grepl, acc), arr.ind=T)
  
  #Check to make sure we don't have more than one match
  if(length(I))
    if(is.integer(I) & length(I)==1) # If only a single acc
      df[1,]=vals[I,]
  else {
    stopifnot(all(!duplicated(I[,'row'])))
    df[I[,'row'],] = vals[I[,'col'], ]
  }
  return(df)
}