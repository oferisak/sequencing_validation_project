
# parse stratification files
stratification_files<-read.table('./data/stratifications.csv',header=F,sep='\t')
colnames(stratification_files)<-c('stratification_name','stratification_file')

# Parse a bed file and output:
# number of regions
# mean length of regions
# total length of regions
get_bed_stats<-function(bed_file){
  bed<-read.table(bed_file,header=F,sep='\t')
  colnames(bed)[1:3]<-c('chr','start','end')
  bed_summary<-bed%>%summarize(n_regions=n(),
                               mean_region_length=mean(end-start),
                               total_region_length=sum(end-start))
  return(bed_summary)
}

strat_files_stats<-NULL
for (i in 1:nrow(stratification_files)){
  strat_name<-stratification_files%>%slice(i)%>%pull(stratification_name)
  strat_file<-stratification_files%>%slice(i)%>%pull(stratification_file)
  info(logger,print(sprintf('Summarizing stratification file for %s: %s',strat_name,strat_file)))
  strat_stats<-get_bed_stats(strat_file)
  strat_files_stats<-strat_files_stats%>%bind_rows(data.frame(stratification_name=strat_name,
                                                              strat_stats,
                                                              stratification_file=strat_file))
}


