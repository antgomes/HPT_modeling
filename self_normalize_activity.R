self_normalize_activity_recipient_label <- function(d_set=q3p11, ind, level="genus",recipient_label=NULL){
  #it normalizes it only in the GEN group (4.1)
  recipient_set = c("B.s.","E.c.","P.a.");
  
  recipient_index=which(recipient_set==recipient_label);
  
  if(length(recipient_index)!=1){
    stop("something went wrong in self_normalize step");
  }
  
  self_strain_set <- c("Bacillus subtilis subtilis 168",
                       "Escherichia coli K-12, MG1655",
                       "Pseudomonas aeruginosa PAO1");
  self_strain = self_strain_set[recipient_index];
  
  ind=d_set$group_field==4.1 & ind;
  key_self_level = d_set[[level]][which(d_set$strain_name==self_strain)[1]];
  
  level_L_activity = d_set$L_activity[d_set[[level]]==key_self_level & ind];
  
  L_activity_mean = mean(level_L_activity);
  L_activity_var = var(level_L_activity);
  
  L_activity_self_normalized=(d_set$L_activity-L_activity_mean)/sqrt(L_activity_var);
  
  return(list(L_activity_self_normalized=L_activity_self_normalized,
              L_activity_mean=L_activity_mean,
              L_activity_var=L_activity_var));
  
}


self_normalize_activity <- function(d_set=q3p11, ind, level="genus",recipient_index=1 ){
  stop("deprecated! Use `self_normalize_activity_recipient_label` instead. Input is `B.s.`, `E.c.`, or `P.a.`.")
  #it normalizes it only in the GEN group (4.1)
  recipient_set = c("B.s.","E.c.","P.a.");
  
  self_strain_set <- c("Bacillus subtilis subtilis 168",
                       "Escherichia coli K-12, MG1655",
                       "Pseudomonas aeruginosa PAO1");
  self_strain = self_strain_set[recipient_index];

  ind=d_set$group_field==4.1 & ind;
  key_self_level = d_set[[level]][which(d_set$strain_name==self_strain)[1]];
  
  level_L_activity = d_set$L_activity[d_set[[level]]==key_self_level & ind];

  L_activity_mean = mean(level_L_activity);
  L_activity_var = var(level_L_activity);
  
  L_activity_self_normalized=(d_set$L_activity-L_activity_mean)/sqrt(L_activity_var);
  
  return(list(L_activity_self_normalized=L_activity_self_normalized,
              L_activity_mean=L_activity_mean,
              L_activity_var=L_activity_var));
  
}

