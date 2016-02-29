fimo_read <- function(fimo_file, sequence_name_field_separator="_-_", ignore_promoter_id_field=F){
  #It returns detailed of 
  rt <- read.csv(fimo_file,header=TRUE, sep="\t");
  
  #For now, I just consider case with `promoter_id=####`. I may have to fix this if I use other fields.
  C <- strsplit(as.character(rt$sequence.name), sequence_name_field_separator);
  n_C <- length(C);
  pre_field <- strsplit(C[[1]],"=");
  field <- pre_field[[1]][1];
  if(ignore_promoter_id_field==F){
    if(field!="promoter_id"){
      stop( "I was expecting the field promoter_id");
    }
  }

  n_fields <- length(C[[1]]);
  for(i in 1:n_fields){
    print(i)
    field_str <- sapply(C, FUN=function(x,p=i){return(x[[p]][1])});
    field_splitted <- strsplit(field_str,"=");
    field_gathered <- unlist(field_splitted);
    field=pre_field[[i]][1];
    if(any(colnames(rt)==field)){
      field=sprintf("%s.1",pre_field[[i]][1]);
    }
    rt[[field]] <- vector(length=n_C); #rt$promoter_id ...
    rt[[field]] <- field_gathered[seq(2,2*n_C,by=2 )];
  }
  if(0){
    stop("Deprecated version! I updated it in Feb/06/2015. This old version was loading only one field");
    rt[[field]] <- vector(length=n_C); #rt$promoter_id ...
    
    #In case there are multiple columns divided by the `sequence_name_field_separator`.
    pre_field1_str <- unlist(C);
    
    #get only `first column` of C.
    field1_str <- pre_field1_str[seq(1,n_fields*n_C,by=n_fields )];
    #Split field1_str by `=`. (promoter_id=#####) 
    field1_splitted <- strsplit(field1_str,"=");
    #Re-gather field1_splitted in promoter_id field
    field1_gathered <- unlist(field1_splitted);
    rt[[field]] <- field1_gathered[seq(2,2*n_C,by=2 )];
  }
  
  return(rt)
  
}

"/home/alg2199/projects/hpt/data/promoter_fasta_set/1_0_1_is_single_TSS=0_motifid=4/fimo.txt"
