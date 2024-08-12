rm(list = ls());  gc()
library(jsonlite)
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(beepr)

v <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
num_temps = 38
Tr = 273.15 + 10
N = 100
Temp_rich = seq(0, num_temps-1, length = num_temps)
k = 0.0000862
temp = Temp_rich +273.15
x = -1/k * (1 /temp - 1/Tr)

# fitted_total = matrix("NA", nrow = 5, ncol = 5)
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
                iter = c(10, 10 ,5 ,5),
                start_lower = c(start_vals_all[1]/100, 0, 0, start_vals_all[4]-10),
                start_upper = c(start_vals_all[1]*100, start_vals_all[2]*2, start_vals_all[3]*2, start_vals_all[4]+10),
                lower = get_lower_lims(long_df$Temp, long_df$Value, model_name = 'sharpeschoolhigh_1981'),
                upper = get_upper_lims(long_df$Temp, long_df$Value, model_name = 'sharpeschoolhigh_1981'),
                supp_errors = 'Y')},
  error = function(e){})

if (!is.null(ss_mod_all)){
  x <- tryCatch(tidy(ss_mod_all), error = function(x){})
  if (!is.null(x)) {
    if (!any(is.nan(x$p.value))) {
      # if (all(x$p.value < 0.05)) {
      output_filename = paste0("../results/",df_name, ".rda")
      saved_data = c(coef(ss_mod_all), AIC(ss_mod_all))
      save(saved_data, file = output_filename)
      # }
    }
  }
}
