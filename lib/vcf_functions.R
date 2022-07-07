# use rtg's vcfstats
rtg_path='/media/SSD/Bioinformatics/Tools/rtg-tools-3.12.1-32d4c2d2/rtg'
get_vcf_stats<-function(vcf_path,target='ALL'){
  message(glue('Getting {vcf_path} stats. Target: {target}'))
  if (target!='ALL'){
    vcfstats_command<-glue('{rtg_path} vcffilter --include-bed={target} -i {vcf_path} -o - | {rtg_path} vcfstats -')
  }else{
    vcfstats_command<-glue('{rtg_path} vcfstats {vcf_path}')
  }
  
  vcfstats_raw_output<-system(vcfstats_command, intern = TRUE)
  vcfstats_parsed<-vcfstats_raw_output%>%stringr::str_split('\\:',simplify = T)%>%trimws()%>%as.data.frame()
  colnames(vcfstats_parsed)<-c('metric','value')
  vcfstats_pivot<-vcfstats_parsed%>%pivot_wider(names_from = metric,values_from = value)
  vcfstats_pivot<-vcfstats_pivot%>%mutate(across(c(`Failed Filters`,
                                                   `Passed Filters`,
                                                   SNPs,
                                                   MNPs,
                                                   Insertions,
                                                   Deletions,
                                                   Indels),
                                                 as.numeric))
  vcfstats_pivot<-vcfstats_pivot%>%mutate(target=target,.before=1)
  return(vcfstats_pivot)
}

get_all_vcf_stats_with_stratifications<-function(){
  all_vcf_stats<-get_samples_sheet_vcf_stats(samples_sheet,target = 'ALL')
  for (i in 1:nrow(stratification_files)){
    strat_name<-stratification_files%>%slice(i)%>%pull(stratification_name)
    message(glue('vcf_functions: Analyzing stratification region: {strat_name}'))
    strat_file<-stratification_files%>%slice(i)%>%pull(stratification_file)
    strat_vcf_stats<-get_samples_sheet_vcf_stats(samples_sheet,strat_file)
    strat_vcf_stats$target<-strat_name
    all_vcf_stats<-all_vcf_stats%>%bind_rows(strat_vcf_stats)
  }
  return(all_vcf_stats)
}

plot_vcfstats_results<-function(all_vcf_stats,target_filter='ALL'){
  toPlot<-all_vcf_stats%>%filter(target==target_filter)%>%select(c(group_name,sample_name,
                                   target,
                                   `Passed Filters`,
                                   SNPs,
                                   Insertions,
                                   Deletions,
                                   Indels))%>%
    pivot_longer(-c(group_name,sample_name,target))%>%
    mutate(name=forcats::fct_relevel(factor(name),'Passed Filters','SNPs','Deletions','Insertions'))
  toPlot$group=stringr::str_replace(toPlot$group_name,'_hg.+','')
  toPlot$ref=stringr::str_match(toPlot$group_name,'hg.+')[,1]
  g<-toPlot%>%ggplot(aes(x=forcats::fct_reorder(factor(sample_name),order(toPlot$group_name)),y=value,fill=group,col=group,label=value))+
    geom_col(width= 0.5,alpha=0.5,position='dodge')+geom_text(vjust=0.5,hjust=-0.5,col='black')+
    facet_grid(name~ref,scales='free')+
    ylim(0,max(toPlot$value)+0.1*max(toPlot$value))+
    coord_flip()+labs(fill=NULL,col=NULL,x=NULL,title=target)+
    theme_minimal()+
    theme(legend.position = 'top')+
    scale_color_nejm()+scale_fill_nejm()
  g
  return(g)
}

get_samples_sheet_vcf_stats<-function(samples_sheet,target=NA){
  samples_sheet_vcf_stats<-NULL
  for (i in 1:nrow(samples_sheet)){
    vcf_path<-samples_sheet%>%slice(i)%>%pull(query_path)
    vcf_name<-samples_sheet%>%slice(i)%>%pull(query_name)
    vcf_group<-samples_sheet%>%slice(i)%>%pull(query_group)
    sample_vcf_stats<-get_vcf_stats(vcf_path,target)
    sample_vcf_stats<-sample_vcf_stats%>%mutate(group_name=vcf_group,
                                                sample_name=vcf_name,
                                                .before=1)
    samples_sheet_vcf_stats<-samples_sheet_vcf_stats%>%bind_rows(sample_vcf_stats)
  }
  return(samples_sheet_vcf_stats)
}


# DEPRACTED
# Read VCF and add hom or het status and transition or transversion annotations
read_vcf<-function(vcf_path){
  vcf_file<-read.vcfR(vcf_path,verbose = F)
  vcf<-vcf_file%>%vcfR2tidy()
  # fix chromosome levels (to be ordered correctly)
  chr_levels<-paste0('chr',c(1:22,'X','Y','M'))
  vcf$fix<-vcf$fix%>%mutate(CHROM=factor(CHROM,levels=chr_levels))
  vcf$fix<-vcf$fix%>%mutate(ti_tv=ifelse((REF!=ALT & nchar(REF)==1 & nchar(ALT)==1),
                                         case_when(REF=='A' & ALT=='G' ~ 'Ti',
                                                   REF=='G' & ALT=='A' ~ 'Ti',
                                                   REF=='C' & ALT=='T' ~ 'Ti',
                                                   REF=='T' & ALT=='C' ~ 'Ti',
                                                   TRUE~'Tv'),NA),
                            het_hom=ifelse(AF=='1.000','hom',ifelse(AF=='0.500','het',NA)))
  return(vcf)
}

# Read all VCFs
read_all_vcfs<-function(samples_sheet){
  all_vcfs<-list()
  for (i in 1:nrow(samples_sheet)){
    vcf_path = samples_sheet%>%slice(i)%>%pull(query_path)
    query_name =samples_sheet%>%slice(i)%>%pull(query_name)
    info(logger,sprintf('Parsing %s VCF..',query_name))
    query_group=samples_sheet%>%slice(i)%>%pull(query_group)
    vcf<-read_vcf(vcf_path)
    if (!(query_group %in% names(all_vcfs))){
      all_vcfs[[query_group]]<-list()
    }
    all_vcfs[[query_group]][[query_name]]<-vcf
  }
  return(all_vcfs)
}

# Parse the VCF fix table and produce met
get_vars_metrics<-function(vcf){
  quantile_probs<-c(0.01,0.05,0.25,0.5,0.75,0.95,0.99)
  quantile_names<-paste0('DP_',100*quantile_probs)
  vcf_metrics <- vcf$fix %>%
    summarize(n = n(),
              DP_mean=mean(DP),
              data.frame(t(quantile(DP, probs = quantile_probs))) %>% `colnames<-`(quantile_names),
              MQ = mean(MQ),
              ti_tv_ratio=sum(ti_tv=='Ti',na.rm = T)/sum(ti_tv=='Tv',na.rm = T),
              het_hom_ratio=sum(het_hom=='het',na.rm = T)/sum(het_hom=='hom',na.rm = T))
  return(vcf_metrics)
}

# produce per-chromosome and all variant stats
analyze_vcf<-function(vcf){
  all_vars_stats_all <- vcf$fix%>%
    get_vars_metrics()%>%mutate(CHROM='All',FILTER='ALL')
  per_chr_stats_all<-vcf$fix%>%group_by(CHROM)%>%
    get_vars_metrics()%>%mutate(FILTER='ALL')
  all_vars_stats_pass <- vcf$fix%>%filter(FILTER=='PASS')%>%
    get_vars_metrics()%>%mutate(CHROM='All',FILTER='PASS')
  per_chr_stats_pass<-vcf$fix%>%filter(FILTER=='PASS')%>%group_by(CHROM)%>%
    get_vars_metrics()%>%mutate(FILTER='PASS')
  overall_var_stats<-all_vars_stats_all%>%bind_rows(all_vars_stats_pass,per_chr_stats_all,per_chr_stats_pass)%>%relocate(CHROM,FILTER)
  return(overall_var_stats)
}

plot_vcfs_coverage<-function(all_vcfs){
  dp_stats<-NULL
  for (query_group in names(all_vcfs)){
    for (query_name in names(all_vcfs[[query_group]])){
      vcf<-all_vcfs[[query_group]][[query_name]]
      # add depth info for all the chromosomes
      dp_stats<-dp_stats%>%bind_rows(data.frame(query_name,query_group,CHROM='ALL',DP=vcf$fix%>%filter(CHROM!='chrM')%>%pull(DP)))
      # add depth info for each chromosome
      dp_stats<-dp_stats%>%bind_rows(data.frame(query_name,query_group,vcf$fix%>%filter(CHROM!='chrM')%>%group_by(CHROM)%>%select(CHROM,DP)))
    }
  }
  coverage_plot<-dp_stats%>%ggplot(aes(y=DP,x=query_name,fill=query_group))+
    geom_violin(alpha=0.5,draw_quantiles = c(0.25,0.5,0.75))+
    #geom_boxplot(outlier.shape = NA,alpha=0.75,width=0.1)+
    theme_minimal()+
    coord_flip()+
    scale_fill_nejm()+
    labs(title='Variant Coverage Per Sample',y='Variant coverage',x=NULL,fill=NULL)+
    theme(legend.position = 'top',panel.grid.minor = element_blank(),plot.title = element_text(face='bold',hjust = 0.5),axis.text=element_text(size=12))+
    scale_y_continuous(breaks=c(10,30,50,100,250))
  return(coverage_plot)
}
