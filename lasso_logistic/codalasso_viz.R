

library(ggplot2)

performance_plot <- function(model){
  
  pred_performance <- data.frame(Lambdas=rep(model$lambdas, 2),
                                        AUC_mean=c(model$train_auc_mean, model$test_auc_mean),
                                        AUC_sd=c(model$train_auc_sd, model$test_auc_sd),
                                        Type=rep(c("Train", "Test"), each=length(model$lambdas)))
  
  pred_performance_plot <- ggplot(pred_performance, aes(x=Lambdas, y=AUC_mean, color=Type)) + 
    geom_point() + geom_line() +
    scale_x_continuous(breaks=model$lambdas)+ylim(0.5, 1)+
    labs(x="Lambda", y="AUC") + 
    theme_bw() + theme(legend.position="top")
  
  var_selection_df <- data.frame(Lambdas=model$lambdas,
                                 Instability = model$instability,
                                 num_features = model$num_features_mean)
  
  instability_plot <- ggplot(var_selection_df, aes(x=Lambdas, y=Instability)) +
    geom_point() + geom_line() + 
    scale_x_continuous(breaks=model$lambdas) +
    labs(x="Lambda", y="Instability") + 
    theme_bw()
  
  num_feature_plot <- ggplot(var_selection_df, aes(x=Lambdas, y=num_features))+
    geom_point() + geom_line() + 
    scale_x_continuous(breaks=model$lambdas) +
    labs(x="Lambda", y="Number of features") + 
    theme_bw()
  
  
  combined_plot <- pred_performance_plot +(instability_plot / num_feature_plot)
  
  return(combined_plot)
  
}



