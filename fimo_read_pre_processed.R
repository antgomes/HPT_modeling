#fimo_read_pre_processed:
#
#It reads fimo file, and pre_process analysis. 
#It outputs `max_L10pvalue` and `n_motif_hits` per promoter_id;
#
#
#

source('/usr/local/home/alg2199/projects/hpt/library/file_parser/fimo_read.R');
source('/usr/local/home/alg2199/projects/hpt/library/Rmysql/get_data_from_hpt_query.R');

fimo_read_pre_processed <- function(fimo_file, motif_thres=3, only_posiive_strand=T, load_q2p1=T){
  #load_q2p1 is TRUE indicates that reference set of`promoter_id` is taken from `q2p1` query.
  
  fr <- fimo_read(fimo_file);
  
  if(only_posiive_strand){
    fr <- fr[fr$strand=="+",];
  }
  fr <- fr[-log10(fr$p.value)>motif_thres,];
  
  if(!load_q2p1){
    stop("error! This script uses q2p1 to match promoter_id oroder;")
  }
  q2p1 <- get_data_from_hpt_query(2.1); #Get set of `promoter_ids`;
  
  
  a_nhits <- aggregate(data.frame(n_hits=-log10(fr$p.value)),
                       by=list(promoter_id=fr$promoter_id),
                       length);
  
  a_maxL10pvalue <- aggregate(data.frame(maxL10pvalue=-log10(fr$p.value)),
                              by=list(promoter_id=fr$promoter_id),
                              max);
  
  m <- merge(q2p1, a_nhits, by="promoter_id",all.x=T);
  m$n_hits[is.na(m$n_hits)]=0;
  m <- merge(m, a_maxL10pvalue, by="promoter_id",all.x=T);
  m$maxL10pvalue[is.na(m$maxL10pvalue)]=2;
  
  return(m);
  
}