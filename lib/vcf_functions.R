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
