# create a dataframe describing the files in the vcfs folder
wd<-getwd()
vcf_files <- list.files('./data/input_vcfs/',pattern = '*vcf$|*.vcf.gz$', recursive = T)
vcf_files <- vcf_files %>% as.data.frame() %>%
  rename('files'='.')%>%
  separate(col=files,into = c('dir','file'),remove=F,sep = '/')%>%
  filter(!grepl('.tbi',file))%>%# filter out vcf indexes
  mutate(files=paste0(sprintf('%s/data/input_vcfs/',wd),files))
