#Aug/11/2015.
#
#I am preparing the figures for the paper: "promoter promiscuity".
#
#This is figure 2:  Describes Biophysics
#         A: Linear regression on   
#           (left)  Activity ~ GC + delta_G + RegulonDB + Sigma70; 
#           (right) Activity ~ Sigma70


#Linear regression on Full model (Activity ~ GC + delta_G + RegulonDB + Sigma70)
#Linear regression on simple model (Activity ~ Sigma70)
#
#

library(MASS);
library(ggplot2);
#library(plyr); #It has the function ddply, which let us compute correlation by `aggregating` per recipient.

source("/usr/local/home/alg2199/projects/hpt/library/Rmysql/get_data_from_hpt_query.R");
source("/usr/local/home/alg2199/projects/hpt/library/antoniostats/closestinterval2.R");
source("/usr/local/home/alg2199/projects/hpt/library/activity_rescaled/self_normalize_activity.R");
source('/usr/local/home/alg2199/projects/hpt/library/antoniostats/cluster_lowcount_bins.R');
source('/usr/local/home/alg2199/projects/hpt/library/file_parser/fimo_read_pre_processed.R');


#Load data

save_path="/usr/local/home/alg2199/projects/hpt/figures/paper_promiscuity/FigureLinearRegress/";
if(!dir.exists(save_path)){
  dir.create(save_path);
}

#query_recipient_set=1;#:3;#This is going to define in which organism the code will be run. I use it for paralel processing.
#recipient_index=1;
recipient_set <- c("B.s.","E.c.", "P.a.");
cur_recipient=recipient_set[recipient_index];
exp_id_set <- c(17, 21, 29); 
#exp_id_set <- c(17, 21, 29, 33, 35);   #33, 35 are Bs and Pa replicate.
total_thres=15;
sTSS_thres=0.8;
do_save=1;
#rna_thres=5;
sTSS_thres=0.8;
motif_thres=3;
rna_thres=5;
dna_thres=1;
n_subsample=500;
#model_integer=1#;+16;#7; #1,2, 4, 7, 15, 23, 36, 45.
#1: GC; 2: delta_G; 4: MaxSigma70; 8=RegulonDB; 16=Genus; #32=n_sigma70_hits


fimo_file_set=c("/usr/local/home/alg2199/projects/hpt/figures/paper_promiscuity/FigureMotif/topn200_mtype1.0_sTSSthres0.80_rnagt5/2014-03-20_LB_B.s./top_wspan-50to0.fimo",
                "/usr/local/home/alg2199/projects/hpt/figures/paper_promiscuity/FigureMotif/topn200_mtype1.0_sTSSthres0.80_rnagt5/2014-03-20_LB_E.c./top_wspan-50to0.fimo",
                "/usr/local/home/alg2199/projects/hpt/figures/paper_promiscuity/FigureMotif/topn200_mtype1.0_sTSSthres0.80_rnagt5/2014-03-20_LB_P.a./top_wspan-50to0.fimo");

recipient_strain_set <- c("Bacillus subtilis subtilis 168",
                          "Escherichia coli K-12, MG1655",
                          "Pseudomonas aeruginosa PAO1");

q2p1 <- get_data_from_hpt_query(2.1);
q2 <- get_data_from_hpt_query(2);

#Loading data.
if(do_load_data){
  q3p12_set <- list();
  q13p01_set <- list();
  fr_set <- list();
  
  for(i in 1:3){
    exp_id = exp_id_set[i];
    
    q3p12_set[[i]] <- get_data_from_hpt_query(3.12,
                                              exp_id);
    ind <- (q3p12_set[[i]]$dna_count+q3p12_set[[i]]$rna_count)>total_thres &
      q3p12_set[[i]]$rna_count>rna_thres &
      q3p12_set[[i]]$TSS_fraction_around_median>sTSS_thres &
      q2p1$is_unique &
      q3p12_set[[i]]$group_field==4.1;
    
    q13p01_set[[i]] <- get_data_from_hpt_query(13.01,
                                               exp_id);
    
    if(any(is.na(q13p01_set[[i]]$delta_G[ind]))){
      stop();
    }
    
    l = self_normalize_activity_recipient_label(q3p12_set[[i]], ind, recipient_label = q3p12_set[[i]]$recipient[1]);
    q3p12_set[[i]]$L_activity_self_normalized=l$L_activity_self_normalized;
    
    fr_set[[i]] <- fimo_read_pre_processed(fimo_file_set[i],
                                           motif_thres = motif_thres,
                                           only_posiive_strand = T);
    
    #Add normalizing part:
  }
}

#Preparing data_set
#     It should contain: L_activtiy_self_normalized, GC_content_bin, delta_G, maxSigma70motif, nSigma70_hits;
if(1){
  for(i in 1:3){
    
    ind_informative <- (q3p12_set[[i]]$rna_count + q3p12_set[[i]]$dna_count)>15 &
      q3p12_set[[i]]$rna_count>rna_thres &
      q3p12_set[[i]]$dna_count>dna_thres &
      q3p12_set[[i]]$group_field==4.1 &
      q3p12_set[[i]]$TSS_fraction_around_median>0.8 & 
      q2p1$is_unique;
    
    d_set_cur <- data.frame(promoter_id=q3p12_set[[i]]$promoter_id,
                            L_activity=q3p12_set[[i]]$L_activity,
                            L_activity_self_normalized=q3p12_set[[i]]$L_activity_self_normalized,
                            GC_content_per_promoter_bin=floor(q3p12_set[[i]]$GC_content_per_promoter/0.1)*0.1,
                            GC_content_per_promoter=q3p12_set[[i]]$GC_content_per_promoter,
                            delta_G_bin=round(q13p01_set[[i]]$delta_G/5)*5,
                            delta_G=q13p01_set[[i]]$delta_G,
                            delta_G_per_bp=(q13p01_set[[i]]$delta_G)/(185-q13p01_set[[i]]$median),
                            maxSigma70motif_bin=floor(fr_set[[i]]$maxL10pvalue/0.5)*0.5,
                            maxSigma70motif=fr_set[[i]]$maxL10pvalue,
                            n_hits=fr_set[[i]]$n_hits,
                            rna_count=q3p12_set[[i]]$rna_count,
                            dna_count=q3p12_set[[i]]$dna_count,
                            genus=q3p12_set[[i]]$genus,
                            recipient=q3p12_set[[i]]$recipient[1]);
    d_set_cur <- d_set_cur[ind_informative,];
    
    if(i==1){
      d_set <- d_set_cur;
    }else{
      d_set <- rbind(d_set,
                     d_set_cur);
    }
  }
}

print("d_set is defined");

#1: GC; 2: delta_G; 4: MaxSigma70; 8=RegulonDB; 16=Genus; #32=n_sigma70_hits
if(model_integer==1){
  model_name="GC_only"; #GC
  
  model <- L_activity ~ (GC_content_per_promoter);
  
}

if(model_integer==2){
  model_name="delta_G_only"; #GC
  
  model <- L_activity ~ (delta_G);
}

if(model_integer==4){
  model_name="sigma_only"; #Sigma70
    
  model <- L_activity ~ (maxSigma70motif);

}

if(model_integer==7){
  
  model_name="full"; #GC + delta_G + Sigma70

  #motif_str <- strsplit(condition_str, "/")[[1]][1];
  model <- L_activity ~ (maxSigma70motif) + 
    GC_content_per_promoter + 
    delta_G;
  
}

if(model_integer==15){
  model_name="full_and_regulonDB"; #GC + delta_G + Sigma70 + RegulonDB;
  pre_motif_number_set= sort( 3:84) + 0.01;  #Sigma factors + TF;
  n_motifs = length(pre_motif_number_set);
  
  motif_str <- paste( sprintf("max_L10p_value%d", floor(pre_motif_number_set) ),
                      collapse=' + ');
  model <- as.formula(paste( "L_activity ~ ",
                      "+",
                      "(maxSigma70motif) + ",
                      "GC_content_per_promoter + ",
                      "delta_G + ",
                      motif_str,
                      sep=" "));
    
  
  #Getting motif information.
  for(j in 1:n_motifs){
    q6p01 <- get_data_from_hpt_query(6.01,exp_id_set[i], pre_motif_number_set[j]);
    q6p01$max_L10p_value[is.na(q6p01$max_L10p_value)]=2;
    
    #q6p01$motif_TSS_distance = q3p12$TSS_median - q6p01$stop; 
    #q6p01$motif_TSS_distance[is.na(q6p01$motif_TSS_distance)]=0.1;
    
    if(any(q3p12_set[[1]]$promoter_id!=q6p01$promoter_id)){
      stop("ERROR! UNEXPECTED UNMATCH!");
    }
    if(j==1){
      df_motif <- data.frame(promoter_id=q6p01$promoter_id,
                               max_L10p_value=q6p01$max_L10p_value);
      #df_motif_TSS_distance <- data.frame(promoter_id=q6p01$promoter_id,
      #                                    motif_TSS_distance=q6p01$motif_TSS_distance);
      df_motif_low_score <- data.frame(promoter_id=q6p01$promoter_id,
                                       low_motif_score=q6p01$max_L10p_value<=4);
    }
    else{
      df_motif <- data.frame(df_motif,
                             max_L10p_value=q6p01$max_L10p_value);
      #df_motif_TSS_distance <- data.frame(df_motif_TSS_distance,
      #                                    motif_TSS_distance=q6p01$motif_TSS_distance);
      df_motif_low_score <- data.frame(df_motif_low_score,
                                       low_motif_score=q6p01$max_L10p_value<=4);
    }
  }
  
  names(df_motif) <- c("promoter_id", 
                       sprintf("max_L10p_value%d",
                               floor(pre_motif_number_set) ));
  
  df_motif[is.na(df_motif)] = 2;
  names(df_motif_low_score) <- c("promoter_id", 
                                 sprintf("low_motif_score%d",
                                         floor(pre_motif_number_set) ));
 # names(df_motif_TSS_distance) <- c("promoter_id", 
  #                                  sprintf("motif_TSS_distance%d",
   #                                         floor(pre_motif_number_set) ));
}

if(model_integer==23){
  model_name="full_and_genus"; #GC + delta_G + Sigma70 + Genus;
  
  model <- L_activity ~ (maxSigma70motif) + 
    GC_content_per_promoter + 
    delta_G +
    genus;
}

if(model_integer==36){
  model_name="MaxSigma_and_nhits"; #GC + delta_G + Sigma70 + Genus;
  
  model <- L_activity ~ (maxSigma70motif) + 
    n_hits;
  
}

if(model_integer==45){
  model_name="full_and_nhits"; #GC + delta_G + Sigma70+n_hits;
  
  model <- L_activity ~ (maxSigma70motif) + 
    GC_content_per_promoter + 
    delta_G +
    n_hits;
}

print("model is loaded");

if(exists("df_motif")){
  d_set$order= seq( length(d_set$promoter_id) );
  d <- merge(d_set, df_motif,by="promoter_id",sort = F);
  d <- d[ order(d$order),];
  if(any(d_set$promoter_id!=d$promoter_id)){
    stop("something is wrong");
  }
} else{
  d <- d_set;
}

ind_informative_core <- d$recipient==cur_recipient;

print("fit");
#Model fit.
fit.model <- lm(model, d[ind_informative_core,]);

#Stepwise Regression:

#step.model_full <- stepAIC(fit.model_full, direction="both");
step.model <- stepAIC(fit.model, direction="both");

#RESULTS:
#THE RESULTS ARE SUMMARIZED IN THE FIELD "ANOVA". 
#step.model_full$anova
step.model$anova
step.model_pred <- predict(step.model,d[ind_informative_core,]);

pre_vs_obs <- data.frame(L_activity_predicted=step.model_pred,
                         L_activity_observed=d$L_activity[ind_informative_core],
                         has_sigma70=as.factor(d$maxSigma70motif[ind_informative_core]>2));

r <- cor(pre_vs_obs$L_activity_predicted, 
         pre_vs_obs$L_activity_observed);

ind_subsample <- sample(1:length(pre_vs_obs$L_activity_predicted), n_subsample);

title_str <- sprintf("%s,%s\ncor=%2.2f\nrna>%d,dna>%d,total_count>%d\nsubsample=%d",
                     cur_recipient,
                     model_name,
                     r,
                     rna_thres, dna_thres, total_thres,
		     n_subsample);

pre_vs_obs$subsample_selected = rep(1,length(pre_vs_obs$L_activity_predicted));
pre_vs_obs$subsample_selected[ind_subsample]=2;
pre_vs_obs$has_sigma70_with_subsample=as.numeric(pre_vs_obs$has_sigma70)-1 + pre_vs_obs$subsample_selected/10; 
pre_vs_obs$has_sigma70_with_subsample=as.factor(pre_vs_obs$has_sigma70_with_subsample);

gg <- ggplot(pre_vs_obs,
             aes(x=L_activity_observed, y=L_activity_predicted, 
                 col=has_sigma70_with_subsample, size=subsample_selected)) +
  geom_point() + 
  ggtitle(title_str) +
  scale_size(range = c(1, 3));

print(gg);

save_str = sprintf("L_activity_%s_predicted_versus_observed_rnagt%d_dnagt%d_totalcountgt%d_nsubsample=%d_%s",
                   model_name,
                   rna_thres,
                   dna_thres,
                   total_thres,
                   n_subsample,
                   cur_recipient);

save_figure <- paste(save_path,
                     save_str,
                     ".eps",
                     sep="");

ggsave(save_figure, gg, width = 8, height = 8);
save_data <- paste(save_path,
                   save_str,
                   ".txt",
                   sep="");
t <- coefficients(step.model);
t_value = data.frame(variable=names(t), value=t);
write.table(t_value,save_data, col.names = T, quote = F, row.names=F);
