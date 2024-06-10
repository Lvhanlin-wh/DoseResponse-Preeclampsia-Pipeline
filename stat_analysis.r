
# this R script is used for major statistical analysis in the manuscript 
# An Analytical Pipeline for Dose-Response effect: Laboratory Tests Assessment and Early Pregnancy Preeclampsia Risk
# authors: XW, HL
# last updated: June 2024

library(dplyr)
library(stringr)
library(ggplot2)
library(mgcv)
library(caret)
library(DescTools)
library(epitools)
library(datawizard)
library(gridExtra)
library(grid)
library(Kendall)
library(trend)
library(nlrr)
library(showtext)
library(rBayesianOptimization)
library(stats)
library(MatchIt)
library(tidyr)
library(scales)
library(segmented)
library(boot)
library(caret)
library(pROC)
library(ROCR)
library(glmnet)
showtext_auto()
options(warn=-1)
options(scipen = 999)

data_path <- "~/Data/"

# auxiliary functions
# calculate second and first derivatives numerically using finite differences
d2_numeric <- function(x,y_prob,gam_model,h) { # h is small increment for finite differences
  dy1 <- predict(gam_model, newdata = data.frame(x = x - h), type = "response")
  dy2 <- predict(gam_model, newdata = data.frame(x = x + h), type = "response")
  d2 <- (dy1 - 2 * y_prob + dy2) / h^2
  return(d2)
}
d1_numeric <- function(x,y_prob,gam_model,h){
  dy1 <- predict(gam_model, newdata = data.frame(x = x + h), type = "response")
  d1 <- (dy1 - y_prob) / h
  return(d1)
}
remove_outliers <- function(x, percentile_left, percentile_right) {
  lower_percentile <- quantile(x, percentile_left)
  upper_percentile <- quantile(x, percentile_right)
  x_filtered <- ifelse(x >= lower_percentile & x <= upper_percentile, x, NA)
  return(x_filtered)
}

# global parameters
percentile_left <- 0.005
percentile_right <- 0.995
spline_type <- 'cr'
grid_step <- 1000  # grid step for pathway A
linear_threshold <- 1.5  # threshold for linearity checking
monotonic_threshold <- 0.6  # threshold for momotonicicity checking tau
p_threshold <- 0.05  # threshold for significance checking
h <- 1e-5  # small increment for finite differences
OR_one_intersection_threshold <- 2  # 2%


# load processed lab data for GZ and GG cohorts 
load(paste(data_path, "GZ_cohort.rdata"))  # all lab records for GZ cohort 
load(paste(data_path, "GZ_lab_list.rdata"))  # 99 labs loaded
load(paste(data_path, "GG_cohort.rdata"))  # all lab records for GG cohort
load(paste(data_path, "GG_lab_list.rdata"))  # 86 labs loaded


# to save all turning points for all labs
P_turning_pt <- data.frame(matrix(ncol = 13, nrow = 0))
names(P_turning_pt) <- c("index","lab","edf_P","tau","p","monotonic","turning_pt_x","turning_pt_y","lower_bdy","upper_bdy","left_gp","right_gp","power")

# to save all risk/protective intervals for all labs
OR_interval <- data.frame(matrix(ncol = 8, nrow = 0))
names(OR_interval) <- c("index","lab","intersection_pt","lower_bdy","upper_bdy","gp_left_size","gp_right_size","power")

                        
# run pipeline for GZ cohort
for (k in 1:nrow(GZ_lab_list)){
  
  current_lab <- GZ_lab_list$lab[k]
  GZ_cohort$current_lab <- (GZ_cohort$lab==eval(current_lab))
  this_lab <- GZ_cohort %>% filter(current_lab==TRUE)
  names(GZ_cohort)[ncol(GZ_cohort)] <- eval(current_lab)
  names(this_lab)[ncol(this_lab)] <- eval(current_lab)
  this_lab$final_lab_num <- remove_outliers(this_lab$lab_num,percentile_left,percentile_right)
  this_lab <- this_lab %>% filter(!is.na(final_lab_num))
  
  if(length(unique(this_lab$final_lab_num))<10){
    print(paste("   There are less than 10 unique values after outlier processing for ", eval(current_lab), ", so exclude it from pipeline", sep=""))
  }else{
    
    #### PATHWAY A - turning point ####
    
    data <- this_lab %>% select(final_lab_num, LABEL)  # LABEL is binary label for PE
    names(data) <- c("x","y")
    gam_model_P <- gam(y ~ s(x, bs = spline_type), data = data, family = binomial())
    edf_P <- pen.edf(gam_model_P)
    # partial effect plot for selected smoother
    plot(gam_model_P, main = paste("Partial Effect Plot for Smoother for Probability of Lab ", eval(current_lab), sep=""),
         xlab = paste("Lab Test for ", eval(current_lab), sep = ""), 
         ylab = paste("s(",eval(current_lab),", edf = ",round(edf_P,2),")",sep=""))
    
    ## prediction to align with the smallest decimal interval 
    min_step <- min(diff(sort(unique(data$x))))
    interval <- (max(data$x)-min(data$x))/min_step
    lab_new <- seq(min(data$x), max(data$x), length.out = min(max(interval, grid_step), grid_step*10))
    
    ## GAM fitting
    GAM_fitting_P <- predict(gam_model_P, data.frame(x = lab_new), type = "response", se.fit = TRUE)
    predicted_P <- data.frame(lab_new = lab_new, y_pred = GAM_fitting_P$fit, lower_ci = GAM_fitting_P$fit - 1.96 * GAM_fitting_P$se.fit, upper_ci = GAM_fitting_P$fit + 1.96 * GAM_fitting_P$se.fit)
    predict_plot_P <- predicted_P %>% ggplot(aes(x = lab_new, y = y_pred)) + geom_line(size=1) +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = 'steelblue') +
      geom_point(aes(x = lab_new, y = lower_ci), color = "royalblue", size = 0.7) +
      geom_point(aes(x = lab_new, y = upper_ci), color = "royalblue", size = 0.7) +
      scale_y_continuous(labels = label_number(accuracy = 0.000001)) +
      labs(x = paste("Lab Test (8-14 Weeks of GA) for", eval(current_lab), sep = " "),
           y = "Probability of PE", title = paste("GAM Model for Lab ", eval(current_lab), " and PE Probability (ind=", k, ")", sep="")) +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 8), axis.title = element_text(size = 6), axis.text = element_text(size = 6))

    ## Mann-Kendall test for monotonicity/trends
    mk <- MannKendall(GAM_fitting_P$fit)  
    sen_slope <- sens.slope(GAM_fitting_P$fit) 
    if(mk$sl[1]<p_threshold & abs(mk$tau[1])>=monotonic_threshold){monotonic_conclusion <- "it shows monotonicity significantly"
    }else{monotonic_conclusion <- "it doesn't show monotonicity significantly"}
    print(paste("   In Mann Kendall test tau=", round(mk$tau[1],2), ", p=", round(mk$sl[1],2), ", sen_slope=", round(as.numeric(sen_slope$estimate),5),", ", monotonic_conclusion, sep = ""))
    
    if(abs(edf_P)>=linear_threshold & mk$sl[1]<p_threshold & abs(mk$tau[1])>=monotonic_threshold){
      
      ## fit (x, first deriv of prob) using GAM
      first_derivatives <- d1_numeric(predicted_P$lab_new, predicted_P$y_pred, gam_model_P, h)
      predicted_P$first_deriv <- first_derivatives
      gam_model_fir <- gam(first_derivatives ~ s(lab_new, bs = spline_type), data = predicted_P)
      edf_fir <- pen.edf(gam_model_fir)
      
      ## prediction of first derivative
      GAM_fitting_deriv1 <- predict(gam_model_fir, data.frame(lab_new = lab_new), type = "response", se.fit = TRUE)
      predicted_deriv1 <- data.frame(lab_new = lab_new, y_pred = GAM_fitting_deriv1$fit,
                                     lower_ci = GAM_fitting_deriv1$fit - 1.96 * GAM_fitting_deriv1$se.fit,
                                     upper_ci = GAM_fitting_deriv1$fit + 1.96 * GAM_fitting_deriv1$se.fit)
      predict_plot_deriv1 <- predicted_deriv1 %>%
        ggplot(aes(x = lab_new, y = y_pred), color = 'black') + geom_line(size=0.5) +
        geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.4, fill = 'darkseagreen1') +
        geom_point(aes(x = lab_new, y = lower_ci), color = "darkgreen", size = 0.3) +
        geom_point(aes(x = lab_new, y = upper_ci), color = "darkgreen", size = 0.3) +
        scale_y_continuous(labels = label_number(accuracy = 0.000001)) +
        labs(#x = paste("Lab Test (8-14 Weeks of GA) for", eval(current_lab), sep = " "),
          x="", y = "First Derivative", title = paste(eval(current_lab), " and First Derivative of PE Probability", sep="")) +
        theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 8), axis.title = element_text(size = 6), axis.text = element_text(size = 6))
      predict_plot_deriv1
      
      ## fit (x, second deriv of prob) using GAM
      second_derivatives <- d2_numeric(predicted_P$lab_new, predicted_P$y_pred, gam_model_P, h)
      predicted_P$second_deriv <- second_derivatives
      gam_model_sec <- gam(second_derivatives ~ s(lab_new, bs = spline_type), data = predicted_P)
      edf_sec <- pen.edf(gam_model_sec)
      
      ## prediction of second derivative
      GAM_fitting_deriv2 <- predict(gam_model_sec, data.frame(lab_new = lab_new), type = "response", se.fit = TRUE)
      predicted_deriv2 <- data.frame(lab_new = lab_new, y_pred = GAM_fitting_deriv2$fit,
                                     lower_ci = GAM_fitting_deriv2$fit - 1.96 * GAM_fitting_deriv2$se.fit,
                                     upper_ci = GAM_fitting_deriv2$fit + 1.96 * GAM_fitting_deriv2$se.fit)
      predict_plot_deriv2 <- predicted_deriv2 %>%
        ggplot(aes(x = lab_new, y = y_pred), color = 'black') + geom_line(size=0.5) +
        geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.4, fill = 'darkseagreen1') +
        geom_point(aes(x = lab_new, y = lower_ci), color = "darkgreen", size = 0.3) +
        geom_point(aes(x = lab_new, y = upper_ci), color = "darkgreen", size = 0.3) +
        scale_y_continuous(labels = label_number(accuracy = 0.000001)) +
        labs(#x = paste("Lab Test (8-14 Weeks of GA) for", eval(current_lab), sep = " "),
          x="",y = "Second Derivative", title = paste(eval(current_lab), " and Second Derivative of PE Probability", sep="")) +
        theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 8), axis.title = element_text(size = 6), axis.text = element_text(size = 6))
      predict_plot_deriv2
      
      ## find uncertainty intervals (UI) associated with turning points
      # find interval associated with max - use upper bound
      max <- max(predicted_deriv2$y_pred) 
      max_ind <- which.max(predicted_deriv2$y_pred) 
      max_pt <- predicted_deriv2$lab_new[max_ind]  
      # find the index and value closest to the max to the left as much as possible, i.e. diff(x) is max and diff(y) is min
      left_elg_ind <- which(predicted_deriv2$lab_new <= max_pt)  
      left_ind <- left_elg_ind[which.min(abs(max - predicted_deriv2$upper_ci[left_elg_ind]))]  
      left_pt <- predicted_deriv2$lab_new[left_ind]  
      # find the index and value closest to the max to the right as much as possible
      right_elg_ind <- which(predicted_deriv2$lab_new >= max_pt)  
      right_ind <- right_elg_ind[which.min(abs(max - predicted_deriv2$upper_ci[right_elg_ind]))] 
      right_pt <- predicted_deriv2$lab_new[right_ind]  
      print(paste("   Estimated max of sec deriv is ",round(max_pt,2), " and associated interval is [",round(left_pt,2),",",round(right_pt,2),"]",sep=""))
      
      # find interval associated with min - use lower bound
      min <- min(predicted_deriv2$y_pred) 
      min_ind <- which.min(predicted_deriv2$y_pred)  
      min_pt <- predicted_deriv2$lab_new[min_ind]  
      # find the index and value closest to the min to the left as much as possible
      left_elg_ind <- which(predicted_deriv2$lab_new <= min_pt)  
      left_ind <- left_elg_ind[which.min(abs(min - predicted_deriv2$lower_ci[left_elg_ind]))]  
      left_pt <- predicted_deriv2$lab_new[left_ind]  
      # find the index and value closest to the min to the right as much as possible
      right_elg_ind <- which(predicted_deriv2$lab_new >= min_pt) 
      right_ind <- right_elg_ind[which.min(abs(min - predicted_deriv2$lower_ci[right_elg_ind]))]  
      right_pt <- predicted_deriv2$lab_new[right_ind]  
      print(paste("   Estimated min of sec deriv is ",round(min_pt,2)," and associated interval is [",round(left_pt,2),",",round(right_pt,2),"]",sep=""))
      
      # max of abs value
      max <- max(abs(predicted_deriv2$y_pred))
      max_ind <- which.max(abs(predicted_deriv2$y_pred))  
      max_pt <- predicted_deriv2$lab_new[max_ind] 
      max_original <- predicted_deriv2$y_pred[max_ind]
      # find the index and value closest to the max to the left as much as possible, i.e. diff(x) is max and diff(y) is min
      left_elg_ind <- which(predicted_deriv2$lab_new <= max_pt) 
      left_ind1 <- left_elg_ind[which.min(abs(max_original - predicted_deriv2$upper_ci[left_elg_ind]))]  
      left_ind2 <- left_elg_ind[which.min(abs(max_original - predicted_deriv2$lower_ci[left_elg_ind]))]  
      left_ind <- min(left_ind1, left_ind2)  
      left_pt <- predicted_deriv2$lab_new[left_ind]  
      # find the index and value closest to the max to the right as much as possible
      right_elg_ind <- which(predicted_deriv2$lab_new >= max_pt) 
      right_ind1 <- right_elg_ind[which.min(abs(max_original - predicted_deriv2$upper_ci[right_elg_ind]))]  
      right_ind2 <- right_elg_ind[which.min(abs(max_original - predicted_deriv2$lower_ci[right_elg_ind]))]  
      right_ind <- max(right_ind1, right_ind2)
      right_pt <- predicted_deriv2$lab_new[right_ind]  
      
      cutoff2 <- round(max_pt,3)
      cutoff_y2 <- round(predicted_P$y_pred[max_ind],4)
      cutoff_lower2 <- round(left_pt,3)
      cutoff_upper2 <- round(right_pt,3)
      print(paste("   Estimated turning point (max of abs of y'') is ",cutoff2," [",cutoff_lower2,",",cutoff_upper2,"]",sep=""))
      
      if(mk$tau[1]>0){
        monotonic <- "increasing"
        subgp_threshold <- 1017
      }else{
        monotonic <- "decreasing"
        subgp_threshold <- 565
      }
      data$group <- ifelse(data$x <= cutoff2, "group1", "group2")
      group_left <- sum(data$group == "group1")
      group_right <- sum(data$group == "group2")
      if(group_left<subgp_threshold | group_right<subgp_threshold){power <- "not satisfied"}else{power <- "satisfied"}
      
      new_row <- c(k,current_lab,round(edf_P,2),round(mk$tau[1],4),round(mk$sl[1],4),monotonic,cutoff2,cutoff_y2,cutoff_lower2,cutoff_upper2,group_left,group_right,power)
      P_turning_pt <- rbind(P_turning_pt, new_row)
      names(P_turning_pt) <- c("index","lab","edf_P","tau","p","monotonic","turning_pt_x","turning_pt_y","lower_bdy","upper_bdy","left_gp","right_gp","power")
      
      turning_interval <- predict_plot_P +
        geom_ribbon(data = subset(predicted_P, lab_new <= cutoff2), aes(ymin=0,ymax=y_pred), fill = '#E6AC00', alpha = 0.5) +
        geom_ribbon(data = subset(predicted_P, lab_new > cutoff2), aes(ymin=0,ymax=y_pred), fill = '#DAB98C', alpha = 0.5) +
        geom_vline(xintercept = cutoff_lower2, linetype = "dotted", color = "black") +
        geom_vline(xintercept = cutoff_upper2, linetype = "dotted", color = "black")

      fig <- grid.arrange(predict_plot_P,
                          predict_plot_deriv1 + geom_vline(xintercept = cutoff2, color = "black"),
                          predict_plot_deriv2 + geom_vline(xintercept = cutoff2, color = "black"), 
                          turning_interval, nrow = 4)
      ggsave(fig, filename = paste("Fig_",eval(current_lab),"_ind",k,".png",sep=""),
             path = paste(result_path, "GZ_pathway_A_",spline_type,"/", sep = ""), device = "png", height = 4, width = 6, units = "in")
    }else{
      print(paste("   Linearity/monotonicity is not satisfied for ", eval(current_lab), ", so exclude it from pipeline", sep=""))
    }

    
    #### PATHWAY B - risk/protective intervals ####
    
    ## calculate dichotomized OR for each x
    threshold_list <- c()
    OR_list <- c()
    for (i in lab_new){  # lab_new is defined in pathway A
      temp <- as.data.frame(table(data$y, data$x>i))
      a <- temp$Freq[4]
      b <- temp$Freq[3]
      c <- temp$Freq[2]
      d <- temp$Freq[1]
      OR <- (a/b)/(c/d)
      threshold_list <- append(threshold_list, i)
      OR_list <- append(OR_list, OR)
    }
    or_df <- data.frame(Column1 = unlist(threshold_list), Column2 = unlist(OR_list))
    names(or_df) <- c("threshold","OR")
    or_df <- or_df %>% filter(!is.na(OR) & OR!=0 & !is.infinite(OR))  # NEW FILTER CONDITION!
    
    ## fit (x, OR) using GAM
    # dataframe "data" is updated, from (lab_value, PE label) to (lab_threshold, OR)
    data <- or_df
    names(data) <- c("x","y")
    gam_model_OR <- gam(y ~ s(x, bs = spline_type), data = data)
    edf_OR <- pen.edf(gam_model_OR)
    # partial effect plot for selected smoother
    plot(gam_model_OR, main = paste("Partial Effect Plot for Smoother for Dichotomized OR of Lab ", eval(current_lab), sep=""),
         xlab = paste("Lab Test for ", eval(current_lab), sep = ""),
         ylab = paste("s(",eval(current_lab),", edf = ",round(edf_OR,2),")",sep=""))
    
    ## prediction of OR
    GAM_fitting_OR <- predict(gam_model_OR, data.frame(x = data$x), type = "response", se.fit = TRUE)
    predicted_OR <- data.frame(threshold_new = data$x, y_pred = GAM_fitting_OR$fit,
                               lower_ci = GAM_fitting_OR$fit - 1.96 * GAM_fitting_OR$se.fit, upper_ci = GAM_fitting_OR$fit + 1.96 * GAM_fitting_OR$se.fit)
    predict_plot_OR <- predicted_OR %>% ggplot(aes(x = threshold_new, y = y_pred)) + geom_line(size=1) +
      geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2, fill = '#008080') +
      geom_point(aes(x = threshold_new, y = lower_ci), color = "#00CED1", size = 0.1) +
      geom_point(aes(x = threshold_new, y = upper_ci), color = "#00CED1", size = 0.1) +
      geom_hline(yintercept=1, linetype="solid", color = "red", size=0.75) +
      labs(x = paste(eval(current_lab), "Lab Test (8-14 Weeks of GA)", sep = " "),
           y = "Odds Ratio of PE", title = paste("GAM Model for Lab ", eval(current_lab), " and PE OR (ind=", k, ")", sep="")) +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 8), axis.title = element_text(size = 6), axis.text = element_text(size = 6))

    ## count number of intersection points with OR=1 and estimation of points
    zero_ind <- abs(predicted_OR$y_pred-1)<=delta  # index of intersection pt, the point closest to OR=1
    zero_pt <- predicted_OR$threshold_new[zero_ind]  # corresponding x for intersection pt
    zero_pt <- unique(round(zero_pt,2))
    
    if(length(zero_pt)!=0 & (100*max(abs(diff(zero_pt)))/(max(predicted_OR$threshold_new)-min(predicted_OR$threshold_new))<=OR_one_intersection_threshold)){
      inter_ind <- (predicted_OR$lower_ci-1<=10*delta & predicted_OR$lower_ci>=1)
      inter_ind_pt <- predicted_OR$threshold_new[inter_ind]
      inter_pt_string_lower <- unique(round(inter_ind_pt,2))
      print(paste("   Lower bounds above OR=1 but very close: ", paste0(inter_pt_string_lower, collapse = ", "), sep = ""))
      inter_ind <- (1-predicted_OR$upper_ci<=delta & predicted_OR$upper_ci<=1)
      inter_ind_pt <- predicted_OR$threshold_new[inter_ind]
      inter_pt_string_upper <- unique(round(inter_ind_pt,2))
      print(paste("   Upper bounds below OR=1 but very close: ", paste0(inter_pt_string_upper, collapse = ", "), sep = ""))
      
      if(length(inter_pt_string_lower)!=0 & length(inter_pt_string_upper)!=0){
        
        min_bdy <- min(inter_pt_string_lower,inter_pt_string_upper)
        max_bdy <- max(inter_pt_string_lower,inter_pt_string_upper)
        
        this_lab$group <- ifelse(this_lab$final_lab_num_win<=min_bdy, "group1", ifelse(this_lab$final_lab_num_win>=max_bdy, "group2", "undefined"))
        group_left <- sum(this_lab$group == "group1")
        group_right <- sum(this_lab$group == "group2")
        undefined <- this_lab %>% filter(group == "undefined") %>% tally()
        print(paste("   Number of lab records that are neither in risk/protective intervals is ",undefined,sep = ""))
        if (group_left<subgp_threshold | group_right<subgp_threshold){power <- "not satisfied"}else{power <- "satisfied"}
        
        new_row <- c(k,current_lab,zero_pt[1],min_bdy,max_bdy,group_left,group_right,power)
        OR_interval <- rbind(OR_interval, new_row)
        names(OR_interval) <-  c("index","lab","intersection_pt","lower_bdy","upper_bdy","gp_left_size","gp_right_size","power")
        
        risk_interval <- predict_plot_OR +
          geom_ribbon(data = subset(predicted_OR, threshold_new <= min_bdy), aes(ymin=0.8*min(predicted_OR$y_pred),ymax=y_pred), fill = '#E6AC00', alpha = 0.5) +
          geom_ribbon(data = subset(predicted_OR, threshold_new > max_bdy), aes(ymin=0.8*min(predicted_OR$y_pred),ymax=y_pred), fill = '#DAB98C', alpha = 0.5)

        ggsave(risk_interval, filename = paste("Fig_",eval(current_lab),"_ind",k,".png",sep=""),
               path = paste(result_path, "GZ_pathway_B_selected/", sep = ""), device = "png", height = 6, width = 5, units = "in")
        
        print(paste("   First interval is [", round(min(predicted_OR$threshold_new),2),",",min_bdy,"]",sep = ""))
        print(paste("   Second interval is [", max_bdy,",", round(max(predicted_OR$threshold_new),2),"]",sep = ""))
        print(paste("   SUMMARY: there is one intersection point with OR=1 for", eval(current_lab)), sep = "")
      }else{
        ggsave(predict_plot_OR, filename = paste("Fig_",eval(current_lab),"_ind",k,".png",sep=""),
               path = paste(result_path, "GZ_pathway_B_non_selected/", sep = ""), device = "png", height = 6, width = 5, units = "in")
        print(paste("   SUMMARY: there are none or more than one intersection point with OR=1 for ", eval(current_lab), ", exclude from pipeline", sep = ""))
      }
    }else{
      ggsave(predict_plot_OR, filename = paste("Fig_",eval(current_lab),"_ind",k,".png",sep=""),
             path = paste(result_path, "GZ_pathway_B_non_selected/", sep = ""), device = "png", height = 6, width = 5, units = "in")
      print(paste("   SUMMARY: there are none or more than one intersection point with OR=1 for ", eval(current_lab), ", exclude from pipeline", sep = ""))
    }
  }
}

P_turning_pt <- P_turning_pt %>% filter(power == "satisfied" & edf_P >= 2)
#OR_interval <- OR_interval %>% filter(power == "satisfied")


# PSM
logit_model <- glm(LABEL ~ age + parity + lab_GA + BMI, data = GZ_cohort, family = binomial(link = "logit"))
propensity_scores <- predict(logit_model, type = "response")
GZ_cohort$propensity_score <- propensity_scores

PSM_delta <- 0.2  # caliper parameter
ratio <- 1  # matching ratio
subgp_threshold <- 565  # relaxed power requirement

matched_ctr_list <- matched_trt_list <- unmatched_ctr_list <- unmatched_trt_list <- c()
RR_est_list <- RR_lower_list <- RR_upper_list <- c()
lab <- trt_PE_list <- trt_nonPE_list <- ctr_PE_list <- ctr_nonPE_list <- c()

for (k in 1:nrow(P_turning_pt)){
  
  current_lab <- P_turning_pt$lab[k]
  cutoff <- as.numeric(P_turning_pt$turning_pt_x[k])
  index <- P_turning_pt$index[k]
  GZ_cohort$current_lab <- (GZ_cohort$lab==eval(current_lab))
  this_lab <- GZ_cohort %>% filter(current_lab==TRUE)
  names(GZ_cohort)[ncol(GZ_cohort)] <- eval(current_lab)
  names(this_lab)[ncol(this_lab)] <- eval(current_lab)
  this_lab$final_lab_num <- remove_outliers(this_lab$lab_num,percentile_left,percentile_right)
  this_lab <- this_lab %>% filter(!is.na(final_lab_num))
  this_lab$lab_flag <- ifelse(this_lab$final_lab_num_win >= cutoff, 1, 0)
  
  matched_data <- matchit(lab_flag ~ propensity_score, data = this_lab, method = 'nearest',  # nearest, optimal
                          ratio = ratio, caliper = PSM_delta*sd(log(propensity_scores)))
  
  if (1 %in% match.data(matched_data)$LABEL && 0 %in% match.data(matched_data)$LABEL){
    
    matched_ctr <- as.data.frame(summary(matched_data)$nn)[4,1]
    matched_trt <- as.data.frame(summary(matched_data)$nn)[4,2]
    unmatched_ctr <- as.data.frame(summary(matched_data)$nn)[5,1]
    unmatched_trt <- as.data.frame(summary(matched_data)$nn)[5,2]
    plot(matched_data, type = "jitter", interactive = FALSE)   # expect not too many unmatched
    
    matched_data <- match.data(matched_data)
    trt_PE <- table(matched_data$LABEL, matched_data$lab_flag)[2,2]
    trt_nonPE <- table(matched_data$LABEL, matched_data$lab_flag)[1,2]
    ctr_PE <- table(matched_data$LABEL, matched_data$lab_flag)[2,1]
    ctr_nonPE <- table(matched_data$LABEL, matched_data$lab_flag)[1,1]
    RR <- riskratio.wald(t(table(matched_data$LABEL, matched_data$lab_flag)))
    RR_est <- as.data.frame(RR$measure)[2,1]  # RR
    RR_lower <- as.data.frame(RR$measure)[2,2]  # lower bound
    RR_upper <- as.data.frame(RR$measure)[2,3]  # upper bound
    
    matched_ctr_list <- append(matched_ctr_list, matched_ctr)
    matched_trt_list <- append(matched_trt_list, matched_trt)
    unmatched_ctr_list <- append(unmatched_ctr_list, unmatched_ctr)
    unmatched_trt_list <- append(unmatched_trt_list, unmatched_trt)
    RR_est_list <- append(RR_est_list, RR_est)
    RR_lower_list <- append(RR_lower_list, RR_lower)
    RR_upper_list <- append(RR_upper_list, RR_upper)
    lab <- append(lab, current_lab)
    trt_PE_list <- append(trt_PE_list, trt_PE)
    trt_nonPE_list <- append(trt_nonPE_list, trt_nonPE)
    ctr_PE_list <- append(ctr_PE_list, ctr_PE)
    ctr_nonPE_list <- append(ctr_nonPE_list, ctr_nonPE)
    
  }else{
    print(paste("  Couldn't match for lab ",current_lab,sep=""))
    matched_ctr_list <- append(matched_ctr_list, NA)
    matched_trt_list <- append(matched_trt_list, NA)
    unmatched_ctr_list <- append(unmatched_ctr_list, NA)
    unmatched_trt_list <- append(unmatched_trt_list, NA)
    RR_est_list <- append(RR_est_list, NA)
    RR_lower_list <- append(RR_lower_list, NA)
    RR_upper_list <- append(RR_upper_list, NA)
    lab <- append(lab, current_lab)
    trt_PE_list <- append(trt_PE_list, NA)
    trt_nonPE_list <- append(trt_nonPE_list, NA)
    ctr_PE_list <- append(ctr_PE_list, NA)
    ctr_nonPE_list <- append(ctr_nonPE_list, NA)
  }
}

P_turning_pt$matched_ctr_P <- matched_ctr_list
P_turning_pt$matched_trt_P <- matched_trt_list
P_turning_pt$unmatched_ctr_P <- unmatched_ctr_list
P_turning_pt$unmatched_trt_P <- unmatched_trt_list
P_turning_pt$trt_PE <- trt_PE_list
P_turning_pt$trt_nonPE <- trt_nonPE_list
P_turning_pt$ctr_PE <- ctr_PE_list
P_turning_pt$ctr_nonPE <- ctr_nonPE_list
P_turning_pt$RR_est_P <- round(RR_est_list,2)
P_turning_pt$RR_lower_P <- round(RR_lower_list,2)
P_turning_pt$RR_upper_P <- round(RR_upper_list,2)
P_turning_pt$RR_CI <- paste(P_turning_pt$RR_est_P," [",P_turning_pt$RR_lower_P,",",P_turning_pt$RR_upper_P,"]",sep="")
P_turning_pt$index <- as.numeric(P_turning_pt$index)
P_turning_pt$turning_pt_x <- as.numeric(P_turning_pt$turning_pt_x)
P_turning_pt <- P_turning_pt %>% filter(!is.na(matched_trt_P))


for (k in 1:nrow(OR_interval)){#nrow(OR_interval)){

  current_lab <- OR_interval$lab[k]
  cutoff <- as.numeric(OR_interval$intersection_pt[k])
  cutoff1 <- as.numeric(OR_interval$lower_bdy[k])
  cutoff2 <- as.numeric(OR_interval$upper_bdy[k])
  index <- as.numeric(OR_interval$index[k])
  GZ_cohort$current_lab <- (GZ_cohort$lab==eval(current_lab))
  this_lab <- GZ_cohort %>% filter(current_lab==TRUE)
  names(GZ_cohort)[ncol(GZ_cohort)] <- eval(current_lab)
  names(this_lab)[ncol(this_lab)] <- eval(current_lab)
  this_lab$final_lab_num <- remove_outliers(this_lab$lab_num,percentile_left,percentile_right)
  this_lab <- this_lab %>% filter(!is.na(final_lab_num))
  this_lab$lab_flag <- ifelse(this_lab$final_lab_num_win >= cutoff2, 0, ifelse(this_lab$final_lab_num_win <= cutoff1, 1, NA))
  table(this_lab$lab_flag, exclude = NULL)
  this_lab %>% filter(is.na(lab_flag)) %>% group_by(LABEL) %>% tally()  # how many PEs are excluded within interval [cutoff1, cutoff2]
  this_lab <- this_lab %>% filter(!is.na(lab_flag))
  
  print(paste("   Sample sizes are ", sum(this_lab$lab_flag == 1), ", ", sum(this_lab$lab_flag == 0), sep=""))
  if(sum(this_lab$lab_flag == 1)<subgp_threshold | sum(this_lab$lab_flag == 0)<subgp_threshold){
    
    print(paste("  (Relaxed) power requirement is not satisfied, exclude",""," from pipeline",sep=""))
    matched_ctr_list <- append(matched_ctr_list, NA)
    matched_trt_list <- append(matched_trt_list, NA)
    unmatched_ctr_list <- append(unmatched_ctr_list, NA)
    unmatched_trt_list <- append(unmatched_trt_list, NA)
    RR_est_list <- append(RR_est_list, NA)
    RR_lower_list <- append(RR_lower_list, NA)
    RR_upper_list <- append(RR_upper_list, NA)
    lab <- append(lab, current_lab)
    trt_PE_list <- append(trt_PE_list, NA)
    trt_nonPE_list <- append(trt_nonPE_list, NA)
    ctr_PE_list <- append(ctr_PE_list, NA)
    ctr_nonPE_list <- append(ctr_nonPE_list, NA)
    
  }else{
    
    matched_data <- matchit(lab_flag ~ propensity_score, data = this_lab, method = 'nearest', # nearest, optimal
                            ratio = ratio, caliper = PSM_delta*sd(log(propensity_scores)))
    
    if (1 %in% match.data(matched_data)$LABEL && 0 %in% match.data(matched_data)$LABEL){
      
      matched_ctr <- as.data.frame(summary(matched_data)$nn)[4,1]
      matched_trt <- as.data.frame(summary(matched_data)$nn)[4,2]
      unmatched_ctr <- as.data.frame(summary(matched_data)$nn)[5,1]
      unmatched_trt <- as.data.frame(summary(matched_data)$nn)[5,2]
      plot(matched_data, type = "jitter", interactive = FALSE)   # expect not too many unmatched
      
      matched_data <- match.data(matched_data)
      trt_PE <- table(matched_data$LABEL, matched_data$lab_flag)[2,2]
      trt_nonPE <- table(matched_data$LABEL, matched_data$lab_flag)[1,2]
      ctr_PE <- table(matched_data$LABEL, matched_data$lab_flag)[2,1]
      ctr_nonPE <- table(matched_data$LABEL, matched_data$lab_flag)[1,1]
      RR <- riskratio.wald(t(table(matched_data$LABEL, matched_data$lab_flag)))
      RR_est <- as.data.frame(RR$measure)[2,1]  # RR
      RR_lower <- as.data.frame(RR$measure)[2,2]  # lower bound
      RR_upper <- as.data.frame(RR$measure)[2,3]  # upper bound
      
      matched_ctr_list <- append(matched_ctr_list, matched_ctr)
      matched_trt_list <- append(matched_trt_list, matched_trt)
      unmatched_ctr_list <- append(unmatched_ctr_list, unmatched_ctr)
      unmatched_trt_list <- append(unmatched_trt_list, unmatched_trt)
      RR_est_list <- append(RR_est_list, RR_est)
      RR_lower_list <- append(RR_lower_list, RR_lower)
      RR_upper_list <- append(RR_upper_list, RR_upper)
      lab <- append(lab, current_lab)
      trt_PE_list <- append(trt_PE_list, trt_PE)
      trt_nonPE_list <- append(trt_nonPE_list, trt_nonPE)
      ctr_PE_list <- append(ctr_PE_list, ctr_PE)
      ctr_nonPE_list <- append(ctr_nonPE_list, ctr_nonPE)
      
    }else{
      
      print(paste("  Couldn't match for lab ",current_lab,sep=""))
      matched_ctr_list <- append(matched_ctr_list, NA)
      matched_trt_list <- append(matched_trt_list, NA)
      unmatched_ctr_list <- append(unmatched_ctr_list, NA)
      unmatched_trt_list <- append(unmatched_trt_list, NA)
      RR_est_list <- append(RR_est_list, NA)
      RR_lower_list <- append(RR_lower_list, NA)
      RR_upper_list <- append(RR_upper_list, NA)
      lab <- append(lab, current_lab)
      trt_PE_list <- append(trt_PE_list, NA)
      trt_nonPE_list <- append(trt_nonPE_list, NA)
      ctr_PE_list <- append(ctr_PE_list, NA)
      ctr_nonPE_list <- append(ctr_nonPE_list, NA)
    }
  }
}

OR_interval$matched_ctr_OR <- matched_ctr_list
OR_interval$matched_trt_OR <- matched_trt_list
OR_interval$unmatched_ctr_OR <- unmatched_ctr_list
OR_interval$unmatched_trt_OR <- unmatched_trt_list
OR_interval$trt_PE <- trt_PE_list
OR_interval$trt_nonPE <- trt_nonPE_list
OR_interval$ctr_PE <- ctr_PE_list
OR_interval$ctr_nonPE <- ctr_nonPE_list
OR_interval$RR_est_OR <- round(RR_est_list,2)
OR_interval$RR_lower_OR <- round(RR_lower_list,2)
OR_interval$RR_upper_OR <- round(RR_upper_list,2)
OR_interval$RR_CI <- paste(OR_interval$RR_est_OR," [",OR_interval$RR_lower_OR,",",OR_interval$RR_upper_OR,"]",sep="")
OR_interval$index <- as.numeric(OR_interval$index)
OR_interval$intersection_pt <- as.numeric(OR_interval$intersection_pt)
OR_interval <- OR_interval %>% filter(!is.na(matched_trt_OR))

OR_interval %>% filter(matched_trt_OR<1017)  # manually check these labs for direction/trends and power requirement
