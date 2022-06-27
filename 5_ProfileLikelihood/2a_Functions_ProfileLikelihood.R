MyPaperTheme <-  theme_set(theme_light(base_size=8))+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        #strip.background = element_rect(fill="grey90", colour="grey90",size=1),
        panel.background = element_rect(colour="black",size=1, fill = NA),
        # strip.text.x = element_text(color='black',margin = margin(.05, 0, .05, 0, "cm")), 
        # strip.text.y = element_text(color='black',margin = margin(.05, 0.05, .05, 0.05, "cm")),
        #axis.line = element_line(colour = 'black'), # This is used to have only the x and y axes of the plot in black.
        legend.text = element_text(size=6, margin = margin(l = 1, unit = "pt")),
        legend.title=element_text(size=8),
        legend.spacing.y = unit(0.01, 'cm'),
        #legend.key.size = unit(0.5, 'cm'),
        legend.key.size = unit(0.5, "lines"),    # key size (unit) # This controls the distance between legend icons
        legend.key.height = NULL,                # key height (unit)
        legend.key.width = NULL,                 # key width (unit)
        #legend.position = c(.95, .95),         # to place legend within the plot
        legend.position="bottom",
        legend.justification = c("centre", "top"),
        legend.box.background = element_rect(fill=NA, colour = NA),
        legend.background=element_blank(), ### ADDED to get rid of grey bg in legend key
        #legend.box.just = "right",
        #legend.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(size = 6, face = "bold"),
        axis.title.x = element_text(size = 8),
        axis.text.y = element_text(size = 6, face = "bold"),
        axis.title.y = element_text(size = 8),
        axis.ticks = element_line(color='black'),
        strip.text.x = element_text(color='black', size=8, face="bold"),#set the size and bold face of panel headings of the facet grid.
        strip.text.y = element_text(color='black', size=8, face="bold"),# Not needed, since we dont have a facet grid panel heading on y axis
        strip.background = element_rect(colour=NA, fill=NA))

Update_param_CI_DF <- function(RDS_file, DF) {
  
  all_params <- c()
  for(i in 1:length(RDS_file)){
    all_params[i] <- RDS_file[[i]]$path
  }
  
  # pathid ## how do I fetch which parameter corresponds to which pathid
  # RDS_file$iMek$residuals[RDS_file$iMek$pathid,]
  for(ii in 1:length(all_params)) {
    param <- all_params[ii]
    
    profile_param <- RDS_file[[param]]$residuals[RDS_file[[param]]$pathid,] #the residual of the model when param is varied
    param_values <- RDS_file[[param]]$explored # the values of the parameter "param" explored. The number of values explored = number of nb_steps
    #plot(param_values,profile_param, type = "l") # residuals Vs values of the parameter. Line plot to check
    
    
    # which index number in this array of residual values corresponds to the minimum residual value
    id <- which(profile_param[] == min(profile_param)) 
    
    # fetch the value of the parameter that lies at this id(from all the values that were explored)
    #RDS_file[[param]]$explored[id]
    
    # from the residuals matrix, fetch the column corresponding to id. 
    #This will give the values of all the 28 parameters at the point where the value of param gave minimum value of residual.
    #RDS_file[[param]]$residuals[,id]
    
    min_res <- min(profile_param)
    
    id <- which(profile_param[] == min(profile_param))
    index <- id[1]
    
    
    # Finding the upper limit of the CI
    RDS_file[[param]]$upper_error_index = NA
    RDS_file[[param]]$upper_pointwise = FALSE
    while(index < length(RDS_file[[param]]$explored)){
      index=index+1
      if(profile_param[index] > min_res+qchisq(0.95,DF)){
        RDS_file[[param]]$upper_error_index=index
        RDS_file[[param]]$upper_pointwise = TRUE
        break
      }
    }
    
    
    RDS_file[[param]]$lower_error_index = NA
    RDS_file[[param]]$lower_pointwise = FALSE
    while(index > 1){
      index=index-1
      if(profile_param[index] > min_res+qchisq(0.95,DF)){ 
        RDS_file[[param]]$lower_error_index=index
        RDS_file[[param]]$lower_pointwise = TRUE
        break
      }
    }
    
  }
  
  return(RDS_file)
  
}

make_pl_plot_paper <- function(aggregate_path_list, my_title = "", xaxis_min=-5, xaxis_max=10) {
  
  
  aggregate_path_list_sorted <- aggregate_path_list %>% 
    group_by(path_name) %>% 
    mutate(max_value = max(abs(value))) %>% 
    arrange(desc(-max_value))
  
  segment_xmin <- min(aggregate_path_list_sorted$lv, na.rm = T)
  segment_xmax <- max(aggregate_path_list_sorted$hv, na.rm = T)
  
  y = c(1.5:((length(aggregate_path_list_sorted$path_name))/2))
  yend = c(1.5:((length(aggregate_path_list_sorted$path_name))/2))
  x = rep(segment_xmin,length(y))
  xend = rep(segment_xmax,length(y))
  
  segment_data = data.frame(y,yend,x,xend)
  
  MyCellLineColours <- case_when(
    grepl("XX",aggregate_path_list_sorted$Cell_type) ~ "#FF0000",
    grepl("XO",aggregate_path_list_sorted$Cell_type) ~ "#0000CD")
  
  path_name_factor <- factor(aggregate_path_list_sorted$path_name, levels=unique(aggregate_path_list_sorted$path_name))
  
  ggplot(aggregate_path_list_sorted,aes(x=value, y=path_name_factor)) +
    geom_vline(xintercept=0,color="grey",lty="dashed", size=0.4)+
    geom_pointrange(aes(x=value, y=path_name_factor, xmin=lv, xmax=hv,colour=Cell_type ),position=ggstance::position_dodgev(height=0.7), fatten = 1, size=0.5) +
    scale_color_manual(values=MyCellLineColours)+
    geom_segment(data = segment_data, aes(x = x, y = y, xend = xend, yend = yend), colour="grey", size= 0.4, alpha = 0.5)+
    ggtitle(my_title)+
    coord_cartesian(xlim=c(xaxis_min,xaxis_max))+ # It is imp to define xlim within coord_cartesian because if done outside, the areas where any data points lie outside of the limits, that entire dataset is removed from the plot. 
    MyPaperTheme+
    theme(panel.grid.minor.y = element_line(colour = "black"),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_text(size=6, face="plain"),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 8, face = "bold"),
          legend.title = element_blank())
  
}

Extract_parameter_CI_paths <- function(Path_Aggregates) {
  
  
  Path_Aggregates_paths_df <- data.frame(Path_Aggregates$paths)
  Path_Aggregates_paths_df$pathname <-  row.names(Path_Aggregates_paths_df)
  
  Path_Aggregates_values_df <- Path_Aggregates_paths_df %>% 
    separate(pathname, into=c("path_name","Cell_type"), sep = " ") %>% 
    mutate(lv=round(lv,3),
           hv=round(hv,3),
           value=round(value,3)) %>% 
    select(Cell_type,path_name,lv,hv,value) %>% 
    mutate(CI= paste0(lv,"  :  ",hv))
  
  return(Path_Aggregates_values_df)
  
}

MakeTable_parameter_CI_paths <- function(Path_Aggregates_values_df) {
  
  Path_Aggregates_values_df_XX <- Path_Aggregates_values_df %>% 
    filter(grepl("XX",Cell_type)) %>% 
    arrange(rev(path_name))
  Path_Aggregates_values_df_XO <- Path_Aggregates_values_df %>% 
    filter(grepl("XO",Cell_type)) %>% 
    arrange(rev(path_name))
  
  Path_Aggregates_values_df_both <- left_join(Path_Aggregates_values_df_XX,Path_Aggregates_values_df_XO, by="path_name", suffix=c(".XX", ".XO"))
  Path_Aggregates_values_df_both <- Path_Aggregates_values_df_both %>% 
    mutate(non_overlap = ifelse((hv.XX<lv.XO|hv.XO<lv.XX),"yes",""))
  
  return(Path_Aggregates_values_df_both)
  
}


