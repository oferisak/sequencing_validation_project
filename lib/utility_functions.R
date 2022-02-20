python_path<-'/usr/bin/python2.7'

validate_happy_input<-function(truth_vcf,query_vcf,reference,fp_bed){
  if (!file.exists(truth_vcf)){stop(glue('validate hap.py input: {truth_vcf} is not found!'))}
  if (!file.exists(query_vcf)){stop(glue('validate hap.py input: {query_vcf} is not found!'))}
  if (!file.exists(reference)){stop(glue('validate hap.py input: {reference} is not found!'))}
  if (!file.exists(fp_bed)){stop(glue('validate hap.py input: {fp_bed} is not found!'))}
}

# create happy command
# vcfeval - whether or not to use the vcfeval engine for comparison, otherwise hap.py uses xcmp
run_happy<-function(happy_path,truth_vcf,query_vcf,reference,fp_bed,output_prefix,
                    target_region=NA,
                    vcfeval=TRUE,
                    vcfeval_path='/media/SSD/Bioinformatics/Tools/rtg-tools-3.12.1-32d4c2d2/rtg',
                    vcfeval_template='/media/SSD/Bioinformatics/Databases/hg19/hg19.SDF',
                    stratification='/media/SSD/Bioinformatics/Projects/sequencing_validation/sequencing_validation_project/data/stratifications.csv',
                    log_file='/media/SSD/Bioinformatics/Projects/sequencing_validation/sequencing_validation_project/logs/project.log'){
  
  validate_happy_input(truth_vcf,query_vcf,reference,fp_bed)
  
  happy_command<-sprintf('%s %s %s "%s" -r %s -f %s -o %s -X -V --no-roc --logfile %s',
                         python_path,
                         happy_path,
                         truth_vcf,
                         query_vcf,
                         reference,
                         fp_bed,
                         output_prefix,
                         log_file)
  if (!is.na(target_region)){
    happy_command<-sprintf('%s -T %s',happy_command,target_region)
    
  }
  if (vcfeval){
    happy_command<-sprintf('%s --engine vcfeval --engine-vcfeval-path %s --engine-vcfeval-template %s',happy_command,vcfeval_path,vcfeval_template)
  }
  if (!is.na(stratification)){
    happy_command<-sprintf('%s --stratification %s',happy_command,stratification)
  }
  
  message(sprintf('Running: %s',happy_command))
  info(logger,happy_command)
  output_dir<-sprintf('./output/%s',output_prefix)
  if (!dir.exists(output_dir)){dir.create(output_dir)}
  setwd(output_dir)
  system(happy_command)
  setwd('../..')
}

# collect happy results
collect_happy_results<-function(samples_sheet){
  happy_outcomes<-NULL
  for (i in 1:nrow(samples_sheet)){
    sample<-samples_sheet%>%slice(i)
    query_group<-sample%>%pull(query_group)
    query_name<-sample%>%pull(query_name)
    happy_output_file<-glue('./output/{query_name}/{query_name}.extended.csv')
    # Parse happy output
    if (!file.exists(happy_output_file)){
      message(glue('{query_name} was not properly analyzed. could not find {happy_output_file}'))
      next
      }
    happy_extended<-read.table(sprintf('./output/%s/%s.extended.csv',query_name,query_name),header=T,sep=',')
    happy_outcomes<-happy_outcomes%>%rbind(data.frame(query_group,query_name,happy_extended))
  }
  return(happy_outcomes)
}

# produce happy summary per group
happy_summary_per_group <- function(happy_results) {
  happy_results %>%
    select(query_group,
           Subset,
           Type,
           Filter,
           Subtype,
           contains('METRIC.')) %>%
    pivot_longer(-c(query_group, Subset, Type, Filter, Subtype), names_to = 'Metric') %>%
    mutate(Metric = str_replace(Metric, 'METRIC.', '')) %>%
    filter(Metric != 'Frac_NA') %>%
    group_by(query_group, Subset, Type, Filter, Subtype, Metric) %>%
    dplyr::summarise(skim_without_charts(value))%>%
    select(-c(skim_type,skim_variable,n_missing,complete_rate))
}

# plot happy results
# facet_by - choose which parameters (in addition to type and metric) to facet by
plot_happy_results <-
  function(happy_results,
           query_group_name=NA,
           subset = c('*'),
           subtype = c('*'),
           filter = c('ALL'),
           var_type = c('SNP','INDEL'),
           facet_by=c()) {
    to_plot<-happy_results %>% 
      select('Type', 'Subtype', 'Filter','Subset', contains(c('query_', 'METRIC.'))) %>%
      filter(Subtype %in% subtype & Filter %in% filter & Subset %in% subset & Type %in% var_type)
    if (!is.na(query_group_name)) {
      to_plot <- to_plot %>% filter(query_group %in% query_group_name)
    }
    
    #facet_formula<-'Type ~ Metric'
    facet_formula<-'.~Metric'
    if (length(facet_by)>0){
      facet_formula<-sprintf('%s + %s',facet_formula,paste0(facet_by,collapse = '+'))
    }
    facet_formula<-as.formula(facet_formula)
    
    g_indel<-to_plot%>%filter(Type=='INDEL')%>%
      pivot_longer(-c(Type, Filter, Subtype, Subset, query_group, query_name),names_to = 'Metric') %>%
      filter(Metric != 'METRIC.Frac_NA') %>%
      mutate(Metric=str_replace(Metric,'METRIC.',''))%>%
      ggplot(aes(x = query_group, y = value,fill=query_group)) +
      geom_boxplot(alpha=0.5) +
      facet_wrap(facet_formula) +
      theme_minimal()+
      ylim(NA,1)+labs(x=NULL,title = ifelse(subtype=='*','Indels',subtype))+
      scale_fill_nejm()+
      coord_flip()+
      guides(fill='none')
    
    if (subtype=='*'){
      g_snp<-to_plot%>%filter(Type=='SNP')%>%
        pivot_longer(-c(Type, Filter, Subtype, Subset, query_group, query_name),names_to = 'Metric') %>%
        filter(Metric != 'METRIC.Frac_NA') %>%
        mutate(Metric=str_replace(Metric,'METRIC.',''))%>%
        ggplot(aes(x = query_group, y = value,fill=query_group)) +
        geom_boxplot(alpha=0.5) +
        facet_wrap(facet_formula) +
        theme_minimal()+
        ylim(NA,1)+labs(x=NULL,title = 'SNP')+
        scale_fill_nejm()+
        coord_flip()+
        guides(fill='none')
        g<-g_snp/g_indel
    }else{
      g<-g_indel
    }
    print(g)
    # g<-to_plot%>%
    #   pivot_longer(-c(Type, Filter, Subtype, Subset, query_group, query_name),names_to = 'Metric') %>%
    #   filter(Metric != 'METRIC.Frac_NA') %>%
    #   mutate(Metric=str_replace(Metric,'METRIC.',''))%>%
    #   ggplot(aes(x = query_group, y = value,fill=query_group)) +
    #   geom_boxplot(alpha=0.5) +
    #   facet_wrap(facet_formula) +
    #   theme_minimal()+
    #   ylim(NA,1)+
    #   scale_fill_nejm()+
    #   coord_flip()+
    #   guides(fill='none')
    # print(g)
  }

# generate sample table 
# get a sample name and the happy results and generate a concise sample summary table
generate_sample_table <- function(happy_results, sample_name) {
  sample_table <- happy_results %>% filter(query_name == sample_name &
                                             Subtype == '*' &
                                             Subset != 'TS_boundary' &
                                             Filter == 'PASS') %>%
    mutate(col_name = sprintf('%s: %s', Subset, Type)) %>%
    arrange(desc(Type)) %>% # sort so that SNPs are before INDELS
    select(-c(Subtype, contains(
      c(
        'query_',
        'TRUTH.',
        'QUERY.UNK',
        'Genotype',
        'QQ',
        'Type',
        'Subset',
        'Filter'
      )
    ))) %>%
    t() %>% as.data.frame()
  colnames(sample_table) <-
    str_replace(sample_table['col_name', ], '\\*', 'ALL')
  sample_table <-
    sample_table[setdiff(rownames(sample_table), 'col_name'), ] %>%
    mutate(across(everything(), as.numeric)) %>%
    mutate(across(where(is.numeric), round, 3))
  return(sample_table)
}

# test if sample was already ran 
# check in the output folder, if there is a folder called 'output_prefix' and if so, if there is a summary file there
is_sample_already_analyzed<-function(output_prefix){
  if (dir.exists(sprintf('./output/%s',output_prefix))){
    if (file.exists(sprintf('./output/%s/%s.summary.csv',output_prefix,output_prefix))){
      info(logger,print(sprintf('%s already analyzed, skipping hap.py stage.',output_prefix)))
      return(TRUE)
    }
  }
  return(FALSE)
}

