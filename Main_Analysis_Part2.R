library(fixest)
library(dplyr)
library(splines)
library(stringr)
library(ggplot2)
library(tibble)
library(lme4)
library(purrr)
### This is for Fig.5 and Supplementary FigS3/TableS3
#Group1a versus Group1b

dat1 <- read.csv('mobility_preparedness.csv')
dat2 <- read.csv('supplies_preparedness.csv')
dat3 <- read.csv('building_preapredness.csv')

preprocess_one <- function(d){
  d %>%
    filter(device > 0, !is.na(visits)) %>%
    mutate(
      GEOID = as.factor(GEOID),
      prep_idx = as.integer(prep_idx),
      prep_w   = as.integer(prep_w),
      coastal_flag = as.integer(coastal_flag),
      ida_ge_33    = as.integer(ida_ge_33),
      experienced_wind_33  = factor(experienced_wind_33, levels = c(1, 2)),
      experienced_rain_150 = factor(experienced_rain_150, levels = c(1, 2))
    ) %>%
    mutate(across(c(max_wind, max_rain),
                  ~ as.numeric(scale(.))))
}

# coastal_flag
subset_coast <- function(d, flag){
  d %>%
    mutate(Coast_group = ifelse(coastal_flag == 1, "Coastal", "Inland")) %>%
    filter(coastal_flag == flag)
}

# ida_ge_33==1 （Ida ≥33 m/s）
subset_ida33 <- function(d){
  d %>% filter(ida_ge_33 == 1)
}

# Fit NGB model
#fit_nb <- function(data){fenegbin(visits ~ prep_idx +prep_idx:max_wind + prep_idx:max_rain + offset(log(device)) | GEOID, cluster = ~ GEOID,data = data)}

#Use the indices as needed (replace line45)
fit_nb <- function(data){fenegbin(visits ~ prep_idx + prep_idx:experienced_wind_33 + prep_idx:experienced_rain_150 + offset(log(device)) | GEOID, cluster = ~ GEOID, data = data)}
#fit_nb <- function(data){fenegbin(visits ~ prep_idx + prep_idx:count_wind_21_33 + prep_idx:count_rain_70_150 + offset(log(device)) | GEOID, cluster = ~ GEOID,data = data)}

#  exp(beta)-1 and 95%CI
extract_exp_ci <- function(model, term_in_model, poi_label, threshold_label, pretty_term){
  ct <- fixest::coeftable(model)
  if(!term_in_model %in% rownames(ct)){
    warning(sprintf("Term '%s' not found in model. Skipped.", term_in_model))
    return(NULL)
  }
  beta <- ct[term_in_model, "Estimate"]
  se   <- ct[term_in_model, "Std. Error"]
  tibble(
    POI       = poi_label,
    Threshold = threshold_label,
    Term      = pretty_term,
    est       = exp(beta)-1,
    lo        = exp(beta - 1.96*se)-1,
    hi        = exp(beta + 1.96*se)-1
  )
}

# Terms need to be consistent with line45
#pull_two_terms <- function(model, poi, thr_label){dplyr::bind_rows(extract_exp_ci(model, "prep_idx:max_wind", poi, thr_label, "Wind"), extract_exp_ci(model, "prep_idx:max_rain", poi, thr_label, "Rain"))}
pull_two_terms <- function(model, poi, thr_label){dplyr::bind_rows(extract_exp_ci(model, "prep_idx:experienced_wind_332", poi, thr_label, "Wind"), extract_exp_ci(model, "prep_idx:experienced_rain_1502", poi, thr_label, "Rain"))}
#pull_two_terms <- function(model, poi, thr_label){dplyr::bind_rows(extract_exp_ci(model, "prep_idx:count_wind_21_33", poi, thr_label, "Wind"), extract_exp_ci(model, "prep_idx:count_rain_70_150", poi, thr_label, "Rain"))}


gas <- preprocess_one(dat1)
gro <- preprocess_one(dat2)
bld <- preprocess_one(dat3)

## All CBGs：Coastal vs Inland =========
# --- Coastal ---
gas_c1 <- subset_coast(gas, 1)
gro_c1 <- subset_coast(gro, 1)
bld_c1 <- subset_coast(bld, 1)

m_gas_c1 <- fit_nb(gas_c1)
m_gro_c1 <- fit_nb(gro_c1)
m_bld_c1 <- fit_nb(bld_c1)

df_c1 <- dplyr::bind_rows(
  pull_two_terms(m_gas_c1, "Gasoline", "All Coastal (flag=1)"),
  pull_two_terms(m_gro_c1, "Grocery",  "All Coastal (flag=1)"),
  pull_two_terms(m_bld_c1, "Building", "All Coastal (flag=1)")
)

# --- Inland ---
gas_c0 <- subset_coast(gas, 0)
gro_c0 <- subset_coast(gro, 0)
bld_c0 <- subset_coast(bld, 0)

m_gas_c0 <- fit_nb(gas_c0)
m_gro_c0 <- fit_nb(gro_c0)
m_bld_c0 <- fit_nb(bld_c0)

df_c0 <- dplyr::bind_rows(
  pull_two_terms(m_gas_c0, "Gasoline", "All Inland (flag=0)"),
  pull_two_terms(m_gro_c0, "Grocery",  "All Inland (flag=0)"),
  pull_two_terms(m_bld_c0, "Building", "All Inland (flag=0)")
)

## Inside Ida ≥33 Group：Coastal vs Inland =========
gas_33 <- subset_ida33(gas)
gro_33 <- subset_ida33(gro)
bld_33 <- subset_ida33(bld)

# --- >33 Coastal ---
gas_33_coast <- subset_coast(gas_33, 1)
gro_33_coast <- subset_coast(gro_33, 1)
bld_33_coast <- subset_coast(bld_33, 1)

m_gas_33_coast <- fit_nb(gas_33_coast)
m_gro_33_coast <- fit_nb(gro_33_coast)
m_bld_33_coast <- fit_nb(bld_33_coast)

df_33_coast <- dplyr::bind_rows(
  pull_two_terms(m_gas_33_coast, "Gasoline", ">33 Coastal"),
  pull_two_terms(m_gro_33_coast, "Grocery",  ">33 Coastal"),
  pull_two_terms(m_bld_33_coast, "Building", ">33 Coastal")
)

# --- >33 Inland ---
gas_33_inland <- subset_coast(gas_33, 0)
gro_33_inland <- subset_coast(gro_33, 0)
bld_33_inland <- subset_coast(bld_33, 0)

m_gas_33_inland <- fit_nb(gas_33_inland)
m_gro_33_inland <- fit_nb(gro_33_inland)
m_bld_33_inland <- fit_nb(bld_33_inland)

df_33_inland <- dplyr::bind_rows(
  pull_two_terms(m_gas_33_inland, "Gasoline", ">33 Inland"),
  pull_two_terms(m_gro_33_inland, "Grocery",  ">33 Inland"),
  pull_two_terms(m_bld_33_inland, "Building", ">33 Inland")
)

n_c1 <- c(dat1 = nrow(gas_c1), dat2 = nrow(gro_c1), dat3 = nrow(bld_c1))
n_c0 <- c(dat1 = nrow(gas_c0), dat2 = nrow(gro_c0), dat3 = nrow(bld_c0))

n_33_coast  <- c(dat1 = nrow(gas_33_coast), dat2 = nrow(gro_33_coast), dat3 = nrow(bld_33_coast))
n_33_inland <- c(dat1 = nrow(gas_33_inland), dat2 = nrow(gro_33_inland), dat3 = nrow(bld_33_inland))

plot_grouped <- function(df, Ns, panel_title){
  d <- df %>%
    dplyr::mutate(
      POI  = factor(POI,  levels = c("Gasoline", "Grocery", "Building"),
                    labels = c("dat1", "dat2", "dat3")),
      Term = factor(Term, levels = c("Wind", "Rain"))
    )
  pd <- position_dodge(width = 0.35)
  
  x_labels <- c(
    dat1 = paste0("N=", scales::comma(Ns["dat1"])),
    dat2 = paste0("N=", scales::comma(Ns["dat2"])),
    dat3 = paste0("N=", scales::comma(Ns["dat3"]))
  )
  
  ggplot(d, aes(x = POI, y = est, color = Term, shape = POI)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(group = Term), position = pd, size = 4) +
    geom_errorbar(aes(ymin = lo, ymax = hi, group = Term),
                  position = pd, width = 0.18) +
    scale_color_manual(values = c(Wind = "red", Rain = "blue"), name = "Effect term") +
    scale_shape_manual(values = c("dat1" = 16, "dat2" = 15, "dat3" = 17), name = "Dataset (POI)") +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_bw(base_size = 18) +
    theme(legend.position = "right", panel.grid.major.x = element_blank())
}

# All CBGs: Coast versus Inland (Fig.S3)
p_c1 <- plot_grouped(df_c1, n_c1, "All samples – Coastal (flag=1)")
p_c0 <- plot_grouped(df_c0, n_c0, "All samples – Inland (flag=0)")

print(p_c1)
print(p_c0)

# Ida ≥33 m/s: Coast versus Inland (Main Fig.5)
p_33_coast  <- plot_grouped(df_33_coast,  n_33_coast,  ">33 m/s only – Coastal")
p_33_inland <- plot_grouped(df_33_inland, n_33_inland, ">33 m/s only – Inland")

print(p_33_coast)
print(p_33_inland)

df_c1
df_c0
df_33_coast
df_33_inland


