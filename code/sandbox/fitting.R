setwd("../../temp_interactions/code")
rm(list = ls());  gc()
library(jsonlite)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(beepr)

num_temps = 38
Tr = 273.15 + 10
N = 100
Temp_rich = seq(0, num_temps-1, length = num_temps)
k = 0.0000862
temp = Temp_rich +273.15
x_v = -1/k * (1 /temp - 1/Tr)

pb <- progress::progress_bar$new(total = 5+20000, # 5*N
                                 clear = FALSE,
                                 format ="[:bar] :percent :elapsedfull")

fitted_total = matrix("NA", nrow = 5, ncol = 5)
for(v in 1:5){
  if(!is.null(pb)){
    pb$tick()
  }
  D_names = names(fromJSON("../data/1com-1.json"))
  dn = D_names[v]
  D = fromJSON("../data/1com-1.json")[[dn]]
  df_name = paste0("fitted_", dn)
  df = as.data.frame(abs(D))
  rm(D);  gc()
  colnames(df) = paste0("S", seq(1,ncol(df),length = ncol(df)))
  df$Temp = temp
  ## fitting all at once
  long_df <- df %>%
    pivot_longer(cols = starts_with("S"),  # Columns to reshape
                 names_to = "α",  # Name of the new column for old column names
                 values_to = "Value")
  start_vals_all <- get_start_vals(long_df$Temp, long_df$Value, model_name = 'sharpeschoolhigh_1981')
  start_vals_all[is.na(start_vals_all)] = 0.0
  start_vals_all[is.infinite(start_vals_all)] = 0.0
  ss_mod_all <- tryCatch({
    nls_multstart(Value ~ sharpeschoolhigh_1981(temp = Temp, r_tref, e, eh, th, tref = 10),
                  data = long_df,
                  iter = c(5, 5 ,4 ,4),
                  start_lower = c(0, 0, 0, start_vals_all[4]-10),
                  start_upper = c(start_vals_all[1]*10, start_vals_all[2]*2, start_vals_all[3]*2, start_vals_all[4]+10),
                  lower = get_lower_lims(long_df$Temp, long_df$Value, model_name = 'sharpeschoolhigh_1981'),
                  upper = get_upper_lims(long_df$Temp, long_df$Value, model_name = 'sharpeschoolhigh_1981'),
                  supp_errors = 'Y')},
    error = function(e){})

  if (!is.null(ss_mod_all)){
    x <- tryCatch(tidy(ss_mod_all), error = function(x){})
    if (!is.null(x)) {
      if (!any(is.nan(x$p.value))) {
        # if (all(x$p.value < 0.05)) {
          fitted_total[v,] <- c(coef(ss_mod_all), AIC(ss_mod_all))
        # }
      }
    }
  }
  rm(long_df);  gc()
  ##
  fitted_matrix = matrix("NA", nrow = (ncol(df)-1), ncol = 5)
  for(i in 1:(ncol(df)-1)){
    cS = paste0("S",i)
    start_vals <- get_start_vals(df$Temp, df[[cS]], model_name = 'sharpeschoolhigh_1981')
    start_vals[is.na(start_vals)] = 0.0
    start_vals[is.infinite(start_vals)] = 0.0
    
    formula_str <- paste0("S", i, " ~ sharpeschoolhigh_1981(temp = Temp, r_tref, e, eh, th, tref = 10)")
    formula <- as.formula(formula_str)
    ss_mod <- tryCatch({
      nls_multstart(formula,
                    data = df,
                    iter = c(5, 5 ,4 ,4),
                    start_lower = c(0, 0, 0, start_vals[4]-10),
                    start_upper = c(start_vals[1]*10, start_vals[2]*2, start_vals[3]*2, start_vals[4]+10),
                    lower = get_lower_lims(df$Temp, df[[cS]], model_name = 'sharpeschoolhigh_1981'),
                    upper = get_upper_lims(df$Temp, df[[cS]], model_name = 'sharpeschoolhigh_1981'),
                    supp_errors = 'Y')},
      error = function(e){})
    
    if (!is.null(ss_mod)){
      x <- tryCatch(tidy(ss_mod), error = function(x){})
      if (!is.null(x)) {
        if (!any(is.nan(x$p.value))) {
          # if (all(x$p.value < 0.05)) {
            fitted_matrix[i,] <- c(coef(ss_mod), AIC(ss_mod))
          # }
        }
      }
    }
  }
  colnames(fitted_matrix) = c("B0", "E", "Ed", "Tp", "AIC")
  write.csv(fitted_matrix, paste0("../data/",df_name, ".csv"), row.names = FALSE)
  rm(fitted_matrix); rm(df); gc()
}
colnames(fitted_total) = c("B0", "E", "Ed", "Tp", "AIC")
write.csv(fitted_total, "../data/fitted_all.csv", row.names = FALSE)


# colnames(fitted_total) = c("B0", "E", "Ed", "Tp", "AIC")
beep(sound = 4, expr = NULL)
