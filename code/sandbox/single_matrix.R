rm(list = ls());  gc()
library(ggplot2)
library(jsonlite)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(beepr)

# v <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

v = 1

num_temps = 38
Tr = 273.15 + 10
N = 100
Temp_rich = seq(0, num_temps-1, length = num_temps)
k = 0.0000862
temp = Temp_rich +273.15
x_v = -1/k * (1 /temp - 1/Tr)

D_names = names(fromJSON("../data/1com-1.json"))
dn = D_names[v]
D = fromJSON("../data/1com-1.json")[[dn]]
df_name = paste0("fitted_", dn)
df = as.data.frame(abs(D))
rm(D);  gc()
colnames(df) = paste0("S", seq(1,ncol(df),length = ncol(df)))
df$Temp = temp

pb <- progress::progress_bar$new(total = ncol(df)-1, # 5*N
                                 clear = FALSE,
                                 format ="[:bar] :percent :elapsedfull")

fitted_matrix = matrix("NA", nrow = (ncol(df)-1), ncol = 5)
for(i in 1:(ncol(df)-1)){
  if(!is.null(pb)){
    pb$tick()
  }
  
  cS = paste0("S",i)
  start_vals <- get_start_vals(df$Temp, df[[cS]], model_name = 'sharpeschoolhigh_1981')
  start_vals[is.na(start_vals)] = 0.0
  start_vals[is.infinite(start_vals)] = 0.0
  start_vals[1] = df[[cS]][11]
  start_vals[2] = max(diff(log(df[[cS]]))/diff(x_v))
  formula_str <- paste0("S", i, " ~ sharpeschoolhigh_1981(temp = Temp, r_tref, e, eh, th, tref = 10)")
  formula <- as.formula(formula_str)
  ss_mod <- tryCatch({
    nls_multstart(S3 ~ sharpeschoolhigh_1981(temp = Temp, r_tref, e, eh, th, tref = 10),
                  data = df,
                  iter = c(10, 10 ,4 ,4),
                  start_lower = c(0, 0, 0, start_vals[4]-10),
                  start_upper = c(start_vals[1]*10, start_vals[2]*2, start_vals[3]*2, start_vals[4]+10),
                  # lower = get_lower_lims(df$Temp, df[[cS]], model_name = 'sharpeschoolhigh_1981'),
                  # upper = get_upper_lims(df$Temp, df[[cS]], model_name = 'sharpeschoolhigh_1981'),
                  supp_errors = 'Y')},
    error = function(e){})
  tidy(ss_mod)
  if (!is.null(ss_mod)){
    x <- tryCatch(tidy(ss_mod), error = function(x){})
    if (!is.null(x)) {
      if (!any(is.nan(x$p.value))) {
        # if (all(x$p.value < 0.05)) {
        fitted_matrix[i,] <- as.numeric(c(coef(ss_mod), AIC(ss_mod)))
        # }
      }
    }
  }
}
colnames(fitted_matrix) = c("B0", "E", "Ed", "Tp", "AIC")
write.csv(fitted_matrix, paste0("../results/",df_name, ".csv"), row.names = FALSE)
beep(sound = 4, expr = NULL)


####### Plotting 
i = 3
cS = paste0("S",i)
long_df <- df %>%
  pivot_longer(cols = starts_with("S"),  # Columns to reshape
               names_to = "ID",  # Name of the new column for old column names
               values_to = "rate")
d <- filter(long_df, ID == cS)

preds <- data.frame(fitted = sharpeschoolhigh_1981(temp = new_temp, r_tref = as.numeric(fitted_matrix[i,1]),
                              e = as.numeric(fitted_matrix[i,2]), eh = as.numeric(fitted_matrix[i,3]),
                              th = as.numeric(fitted_matrix[i,1]), tref = 10),
                    temp = new_temp)

# Plot using ggplot2
ggplot(d, aes(temp, rate)) +
  geom_point() +
  geom_line(aes(temp, fitted), preds, col = 'blue') +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = '|α|',
       title = '')
