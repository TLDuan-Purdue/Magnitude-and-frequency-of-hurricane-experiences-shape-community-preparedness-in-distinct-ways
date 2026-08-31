library(fixest)
library(dplyr)
library(splines)
library(stringr)
library(ggplot2)
library(tibble)
library(lme4)
library(tibble)
library(purrr)

# Main results for:
#Associations between hurricane experience and preparedness for mobility needs, daily supplies, and structural reinforcement: 
#All CBGs (Main.3); CBGs with Ida forecast wind greater (and smaller) than hurricane force (33m/s): Group1 versus Group2 (Main Fig.4)

dat1 <- read.csv('mobility_preparedness.csv')
dat2 <- read.csv('supplies_preparedness.csv')
dat3 <- read.csv('building_preapredness.csv')

## ========= 0) Data pre-processing
preprocess_one <- function(d){
  d %>%
    filter(device > 0, !is.na(visits)) %>%
    mutate(
      GEOID = as.factor(GEOID),
      prep_idx = as.integer(prep_idx), prep_w = as.integer(prep_w),ida_ge_33 = as.integer(ida_ge_33),
      experienced_wind_33 = factor(experienced_wind_33, levels = c(1, 2)),
      experienced_rain_150 = factor(experienced_rain_150, levels = c(1, 2))
    ) %>% mutate(across(c(max_wind, max_rain), ~ as.numeric(scale(.))))
}

subset_ida <- function(d, thr){ d %>% mutate(IDA_group = ifelse(Ida0828 > thr, "Ida_1", "Ida_0")) %>% filter(IDA_group == "Ida_1") }

fit_nb <- function(data){fenegbin(visits ~ prep_idx + prep_idx:max_wind + prep_idx:max_rain + offset(log(device)) | GEOID, cluster = ~ GEOID,data = data)}

#Use the indices as needed (replace line33)
#fit_nb <- function(data){fenegbin(visits ~ prep_idx + prep_idx:experienced_wind_33 + prep_idx:experienced_rain_150 + offset(log(device)) | GEOID, cluster = ~ GEOID, data = data)}
#fit_nb <- function(data){fenegbin(visits ~ prep_idx + prep_idx:count_wind_21_33 + prep_idx:count_rain_70_150 + offset(log(device)) | GEOID, cluster = ~ GEOID,data = data)}

# —— Function：Select a group of indice，calculate coef and 95%CI —— 
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
    Term      = pretty_term,                  # e.g., "Wind" / "Rain"
    est       = exp(beta)-1,
    lo        = exp(beta - 1.96*se)-1,
    hi        = exp(beta + 1.96*se)-1
  )
}

# Terms need to be consistent with line33
pull_two_terms <- function(model, poi, thr_label){dplyr::bind_rows(extract_exp_ci(model, "prep_idx:max_wind", poi, thr_label, "Wind"), extract_exp_ci(model, "prep_idx:max_rain", poi, thr_label, "Rain"))}
#pull_two_terms <- function(model, poi, thr_label){dplyr::bind_rows(extract_exp_ci(model, "prep_idx:experienced_wind_332", poi, thr_label, "Wind"), extract_exp_ci(model, "prep_idx:experienced_rain_1502", poi, thr_label, "Rain"))}
#pull_two_terms <- function(model, poi, thr_label){dplyr::bind_rows(extract_exp_ci(model, "prep_idx:count_wind_21_33", poi, thr_label, "Wind"), extract_exp_ci(model, "prep_idx:count_rain_70_150", poi, thr_label, "Rain"))}


## ========= 1) Three types of POIs
gas <- preprocess_one(dat1)        # Gasoline
gro <- preprocess_one(dat2)        # Grocery
bld <- preprocess_one(dat3)        # Building

## ========= 2) Three sampling strategies: All / >=33 / <33 =========

# All CBGs
m_gas_all <- fit_nb(gas)
m_gro_all <- fit_nb(gro)
m_bld_all <- fit_nb(bld)

df_all <- dplyr::bind_rows(
  pull_two_terms(m_gas_all, "Gasoline", "All"),
  pull_two_terms(m_gro_all, "Grocery",  "All"),
  pull_two_terms(m_bld_all, "Building", "All")
)

# =========================
# Group 1: Ida wind >= 33 m/s
# =========================
gas_33_high <- gas %>% filter(ida_ge_33 == 1)
gro_33_high <- gro %>% filter(ida_ge_33 == 1)
bld_33_high <- bld %>% filter(ida_ge_33 == 1)

m_gas_33_high <- fit_nb(gas_33_high)
m_gro_33_high <- fit_nb(gro_33_high)
m_bld_33_high <- fit_nb(bld_33_high)

df_33_high <- dplyr::bind_rows(
  pull_two_terms(m_gas_33_high, "Gasoline", ">=33 m/s"),
  pull_two_terms(m_gro_33_high, "Grocery",  ">=33 m/s"),
  pull_two_terms(m_bld_33_high, "Building", ">=33 m/s")
)

# =========================
# Group 2: Ida wind < 33 m/s
# =========================
gas_33_low <- gas %>% filter(ida_ge_33 == 0)
gro_33_low <- gro %>% filter(ida_ge_33 == 0)
bld_33_low <- bld %>% filter(ida_ge_33 == 0)

m_gas_33_low <- fit_nb(gas_33_low)
m_gro_33_low <- fit_nb(gro_33_low)
m_bld_33_low <- fit_nb(bld_33_low)

df_33_low <- dplyr::bind_rows(
  pull_two_terms(m_gas_33_low, "Gasoline", "<33 m/s"),
  pull_two_terms(m_gro_33_low, "Grocery",  "<33 m/s"),
  pull_two_terms(m_bld_33_low, "Building", "<33 m/s")
)

# Sample sizes
n_all <- c(dat1 = nrow(gas),dat2 = nrow(gro),dat3 = nrow(bld))
n_33_high <- c(dat1 = nrow(gas_33_high),dat2 = nrow(gro_33_high),dat3 = nrow(bld_33_high))
n_33_low <- c(dat1 = nrow(gas_33_low),dat2 = nrow(gro_33_low),dat3 = nrow(bld_33_low))

#Plot function
plot_grouped <- function(df, Ns, panel_title){
  d <- df %>%
    dplyr::mutate(
      POI  = factor(POI,  levels = c("Gasoline", "Grocery", "Building"),
                    labels = c("dat1", "dat2", "dat3")),
      Term = factor(Term, levels = c("Wind", "Rain"))
    )
  pd <- position_dodge(width = 0.35)
  
  # Labels
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
    scale_y_continuous(
      breaks = function(lims) {
        lo <- floor(lims[1] * 100)
        hi <- ceiling(lims[2] * 100)
        seq(lo, hi, by = 2) / 100   # 
      },
      labels = function(x) paste0(round(x * 100), "%")
    )+
    labs(x = NULL, y = NULL, title = NULL) +
    theme_bw(base_size = 18) +
    theme(legend.position = "right", panel.grid.major.x = element_blank())
}

# Plots
# Plots
p_all <- plot_grouped(
  df_all,
  n_all,
  "All samples"
)

p_33_high <- plot_grouped(
  df_33_high,
  n_33_high,
  "Ida forecast wind >= 33 m/s"
)

p_33_low <- plot_grouped(
  df_33_low,
  n_33_low,
  "Ida forecast wind < 33 m/s"
)

print(p_all)
print(p_33_high)
print(p_33_low)

# Output estimates
df_all
df_33_high
df_33_low

