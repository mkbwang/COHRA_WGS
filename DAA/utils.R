DAA <- function(title, 
                source=c("saliva", "plaque"),
                timestamp=c("yr1", "yr2"),
                type=c("taxa", "ko"),
                covariate=c("inc25", "inc75", "Cigarettes", "region", "higherED"),
                prevalence_filter=0.1,
                raw_pval=0.05,
                adj_pval=0.1){
  
  metadata <- read.table(sprintf("counts_cleaning/strata/metadata_%s_%s.tsv", source, timestamp),
                         header=T, sep='\t')
  
  counts <- read.table(sprintf("counts_cleaning/strata/%s_%s_counts_%s.tsv", source, type, timestamp), 
                       header=T,
                       sep='\t')
  
  output <- maaslin3(input_data=counts,
                     input_metadata=metadata,
                     output=sprintf("DAA/%s/%s/%s/%s",timestamp,source,type,covariate),
                     formula = sprintf("~%s",covariate),
                     min_prevalence=prevalence_filter,
                     subtract_median=T,
                     evaluate_only="abundance",
                     warn_prevalence = F,
                     plot_summary_plot = F)
  
  fit_result <- output$fit_data_abundance$results
  fit_result_significant <- fit_result %>% filter(pval_individual < raw_pval) %>%
    select(feature, coef, qval_individual, pval_individual)
  fit_result_significant$significant <- fit_result_significant$qval_individual < adj_pval
  
  # prepare ggplot object
  fit_result_significant$feature <- factor(fit_result_significant$feature,
                                           levels=rev(fit_result_significant$feature))
  
  if (type == "taxa"){
    plot_title <- sprintf("DA taxa of %s samples collected at %s", source, timestamp)
  } else{
    plot_title <- sprintf("DA KEGG Ortholog of %s samples collected at %s", source, timestamp)
  }
  
  
  if(length(unique(fit_result_significant)) == 1){ # no signals after p value adjustment
    visualization <- ggplot(fit_result_significant, aes(x=feature, y=coef, color="#111111")) + 
      geom_bar(stat = "identity", fill = "gray") + 
      xlab("") + ylab(sprintf("Log 2 Abundance Fold Change (%s)", title))
    theme_bw() + coord_flip()+
      ggtitle(plot_title)+
      geom_vline(xintercept = 0, linetype="dashed")
  } else{ # with signals after p value adjustment
    visualization <- ggplot(fit_result_significant, aes(x=feature, y=coef, color=significant)) + 
      geom_bar(stat = "identity", fill = "gray") + 
      xlab("") + ylab(sprintf("Log 2 Abundance Fold Change (%s)", title))+
      scale_color_manual(values=c("#111111" , "#ff3300")) +
      theme_bw() + coord_flip()+guides(color="none")+
      ggtitle(plot_title)+
      geom_vline(xintercept = 0, linetype="dashed")
  }
  
  return(list(signals=fit_result_significant, visual=visualization))
  
}
