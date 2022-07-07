validate_input_bed_file<-function(input_bed_file){
  input_bed<-readr::read_delim(input_bed_file,col_names = c('chr','start','end'))
  if (ncol(input_bed)>3){message(glue('There are more than 3 columns in input bed file {input_bed_file}, will remove them'))}
  
  # if the input was gzipped, save it unpacked then gzip it
  if (grepl('.gz',input_bed_file)){
    #write.table(as.data.frame(input_bed[,1:3]),file=stringr::str_replace(input_bed_file,'.gz',''),row.names = F,col.names = F,sep='\t',quote = F)
    write.table(format(as.data.frame(input_bed[,1:3]),scientific=F,trim=T,justify='none'),file=stringr::str_replace(input_bed_file,'.gz',''),row.names = F,col.names = F,sep='\t',quote = F)
    system(glue('gzip -f {stringr::str_replace(input_bed_file,".gz","")}'))
  }else{
    write.table(format(as.data.frame(input_bed[,1:3]),scientific=F,trim=T,justify='none'),file=stringr::str_replace(input_bed_file,'.gz',''),row.names = F,col.names = F,sep='\t',quote = F)
  }
  
}

# parse stratification files
stratification_files<-list.files('/media/SSD/Bioinformatics/Projects/sequencing_validation/sequencing_validation_project/data/stratifications',full.names = T,recursive = T)
for (strat_file in stratification_files){
  validate_input_bed_file(strat_file)
}

stratification_files<-data.frame(stratification_file=stratification_files,stratification_name=stringr::str_replace(basename(stratification_files),'.bed',''))
# generate the stratifications file
write.table(stratification_files%>%select(stratification_name,stratification_file),
            sep = '\t',
            file='/media/SSD/Bioinformatics/Projects/sequencing_validation/sequencing_validation_project/data/stratifications.csv',
            row.names = F,col.names = F,quote = F)

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


