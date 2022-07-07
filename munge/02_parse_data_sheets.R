# Parse samples sheet ####
#samples_sheet<-read.table('./data/samples_sheet.csv',header=T,sep=',',stringsAsFactors = F)

# use the vcfs to update full file path
# message('Updating samples sheet with full vcf file paths..')
# query_path<-c()
# for (i in 1:nrow(samples_sheet)){
#   sample<-samples_sheet%>%slice(i)
#   query_file<-sample%>%pull(query_file)
#   query_group<-sample%>%pull(query_group)
#   if (!(query_file %in% vcf_files$file)){stop(sprintf('There is no input vcf file with the name %s',query_file))}
#   sample_path<-vcf_files%>%filter(dir==query_group & file==query_file)%>%pull(files)
#   query_path=c(query_path,sample_path)
# }
# samples_sheet$query_path<-query_path

# Parse GIAB sheet ####
giab_sheet<-read.table('./data/giab_sheet.csv',header=T,sep='\t',stringsAsFactors = F)
