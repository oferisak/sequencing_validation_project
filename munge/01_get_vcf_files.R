# create a dataframe describing the files in the vcfs folder
wd<-getwd()
vcf_files <- list.files('./data/input_vcfs/',pattern = '*vcf$|*.vcf.gz$', recursive = T)
samples_sheet <- vcf_files %>% as.data.frame() %>%
  rename('files'='.')%>%
  separate(col=files,into = c('query_group','query_file'),remove=T,sep = '/')%>%
  filter(!grepl('.tbi',query_file))%>%# filter out vcf indexes
  mutate(query_path=glue('{wd}/data/input_vcfs/{query_group}/{query_file}'),
         query_name=stringr::str_extract(query_file,'[^\\.]+'),
         truth=case_when(
           grepl('hg002|21-063|063FILES|21_063',query_name) ~ 'hg002',
           grepl('hg003|21-064|064FILES|21_064',query_name) ~ 'hg003',
           grepl('hg004|21-065|065FILES|21_065',query_name) ~ 'hg004',
         ))
  if (sum(is.na(samples_sheet$truth))>0){stop('Could not identify the relevant GIAB reference for some of the input vcf files')}

message(paste0(c(glue('Collected {nrow(samples_sheet)} vcfs from the ./data/input_vcfs folder'),
               glue('There are {length(unique(samples_sheet$query_group))} groups in the folder')),
               collapse='\n'))
