library(fixest)
library(dplyr)
library(ggplot2)
library(tibble)
library(purrr)

# This is for the main results about property damage (Fig6)
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
    ) 
}



fit_nb <- function(data){
  fenegbin(
    visits ~ prep_idx + prep_idx:PropertyDmgPerCapita + offset(log(device)) | FIPS,
    cluster = ~ FIPS,
    data = data
  )
}

extract_exp_ci <- function(model, term_in_model, poi_label, group_label){
  ct <- fixest::coeftable(model)
  if(!term_in_model %in% rownames(ct)){
    warning(sprintf("Term '%s' not found in model. Skipped.", term_in_model))
    return(NULL)
  }
  beta <- ct[term_in_model, "Estimate"]
  se   <- ct[term_in_model, "Std. Error"]
  tibble(
    POI   = poi_label,
    Group = group_label,
    est   = exp(beta) - 1,
    lo    = exp(beta - 1.96 * se) - 1,
    hi    = exp(beta + 1.96 * se) - 1
  ) 
}

pull_one_term <- function(model, poi, group_label){
  extract_exp_ci(model, "prep_idx:PropertyDmgPerCapita", poi, group_label)
}


gas <- preprocess_one(dat1)
gro <- preprocess_one(dat2)
bld <- preprocess_one(dat3)

# Five different samplings (different hazard scenarios)
get_scenario_data <- function(flag_id){
  if(flag_id == "all"){
    list(gas = gas, gro = gro, bld = bld,
         label = "All samples")
  } else if(flag_id == "ida_gt33"){
    list(gas = gas %>% filter(ida_ge_33 == 1),
         gro = gro %>% filter(ida_ge_33 == 1),
         bld = bld %>% filter(ida_ge_33 == 1),
         label = "Ida > 33 m/s")
  } else if(flag_id == "ida_lt33"){
    list(gas = gas %>% filter(ida_ge_33 == 0),
         gro = gro %>% filter(ida_ge_33 == 0),
         bld = bld %>% filter(ida_ge_33 == 0),
         label = "Ida < 33 m/s")
  } else if(flag_id == "ida_gt33_coast"){
    list(gas = gas %>% filter(ida_ge_33 == 1, coastal_flag == 1),
         gro = gro %>% filter(ida_ge_33 == 1, coastal_flag == 1),
         bld = bld %>% filter(ida_ge_33 == 1, coastal_flag == 1),
         label = "Ida > 33 m/s (Coastal)")
  } else if(flag_id == "ida_gt33_inland"){
    list(gas = gas %>% filter(ida_ge_33 == 1, coastal_flag == 0),
         gro = gro %>% filter(ida_ge_33 == 1, coastal_flag == 0),
         bld = bld %>% filter(ida_ge_33 == 1, coastal_flag == 0),
         label = "Ida > 33 m/s (Inland)")
  } else {
    stop("Unknown scenario id.")
  }
}

scenario_ids <- c("all", "ida_gt33", "ida_lt33", "ida_gt33_coast", "ida_gt33_inland")

#
plot_grouped <- function(df, Ns, panel_title){
  
  d <- df %>%
    mutate(
      POI = factor(POI, levels = c("Gasoline", "Grocery", "Building"),
                   labels = c("dat1", "dat2", "dat3"))
    )
  
  x_labels <- c(
    dat1 = paste0("N=", scales::comma(Ns["dat1"])),
    dat2 = paste0("N=", scales::comma(Ns["dat2"])),
    dat3 = paste0("N=", scales::comma(Ns["dat3"]))
  )
  
  ggplot(d, aes(x = POI, y = est, shape = POI)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(color = "black", size = 4) +
    geom_errorbar(aes(ymin = lo, ymax = hi),
                  width = 0.18, color = "black") +
    scale_shape_manual(
      values = c("dat1" = 16, "dat2" = 15, "dat3" = 17),
      name = "Dataset (POI)"
    ) +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    # scale_y_continuous(
    #   breaks = function(lims) {
    #     lo <- floor(lims[1] * 100)
    #     hi <- ceiling(lims[2] * 100)
    #     seq(lo, hi, by = 3) / 100   
    #   },
    #   labels = function(x) paste0(round(x * 100), "%")
    # )+
    labs(x = NULL, y = NULL, title = NULL) +
    theme_bw(base_size = 18)
    #theme(
    #legend.position = "right",
    #panel.grid.major.x = element_blank()
    #)
}

## ========= Design loops
plots <- list()
dfs   <- list()

for(sid in scenario_ids){
  
  sc <- get_scenario_data(sid)
  
  m_gas <- fit_nb(sc$gas)
  m_gro <- fit_nb(sc$gro)
  m_bld <- fit_nb(sc$bld)
  
  df_sc <- bind_rows(
    pull_one_term(m_gas, "Gasoline", sc$label),
    pull_one_term(m_gro, "Grocery",  sc$label),
    pull_one_term(m_bld, "Building", sc$label)
  )
  
  # N
  Ns <- c(dat1 = nrow(sc$gas), dat2 = nrow(sc$gro), dat3 = nrow(sc$bld))
  
  # Plots
  p_sc <- plot_grouped(df_sc, Ns, sc$label)
  
  plots[[sid]] <- p_sc
  dfs[[sid]]   <- df_sc
}

## Generate Graphs
print(plots[["all"]])
print(plots[["ida_gt33"]])
print(plots[["ida_lt33"]])
print(plots[["ida_gt33_coast"]])
print(plots[["ida_gt33_inland"]])

## Generate Images
dfs[["all"]]
dfs[["ida_gt33"]]
dfs[["ida_lt33"]]
dfs[["ida_gt33_coast"]]
dfs[["ida_gt33_inland"]]


