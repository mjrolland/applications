# Plot functions
# M. Rolland
# 28/05/202

# List of compounds and their label
comp_labels <- function(){
  labs <- c(
    "MEPA_total" = "Metylparab\u00e8ne", 
    "ETPA_total" = "Ethylparab\u00e8ne", 
    "PRPA_total" = "Propylparab\u00e8ne", 
    "BUPA_total" = "Butylparab\u00e8ne", 
    "BPA_total"  = "Bisph\u00e9nol A", 
    "BPS_total"  = "Bisph\u00e9nol S", 
    "BPF_total"  = "Bisph\u00e9nol F", 
    "BPB_total"  = "Bisph\u00e9nol B", 
    "BPAF_total" = "Bisph\u00e9nol AF", 
    "OXBE_total" = "Benzoph\u00e9none 3", 
    "TRCS_total" = "Triclosan", 
    "TRCB_total" = "Triclocarban",
    "MEP"        = "MEP",
    "MnBP"       = "MnBP",
    "MiBP"       = "MiBP",
    "MBzP"       = "MBzP",
    "MEHP"       = "MEHP",
    "MEHHP"      = "MEHHP",
    "MEOHP"      = "MEOHP",
    "MECPP"      = "MECPP",
    "MMCHP"      = "MMCHP",
    "ohMPHP"     = "oh-MPHP",
    "ohMiNP"     = "oh-MiNP",
    "oxoMiNP"    = "oxo-MiNP",
    "cxMiNP"     = "cx-MiNP",
    "ohMINCH"    = "oh-MINCH",
    "oxoMINCH"   = "oxo-MINCH",
    "DEHPms2"    = "DEHP",
    "DiNPms2"    = "DiNP",
    "DINCHms2"   = "DINCH"
  )
  return(labs)
}

# function to plot cohort distribution + individual levels for a given set of
# compounds
plot_hist <- function(samples, comp_group, id, my_t2, my_t3){
  # prepare whole population data
  plot_data <- samples %>%
    filter(compound %in% comp_group) %>%
    group_by(compound) %>%
    mutate(
      q1  = quantile(val, 0.25, na.rm = TRUE),
      q3  = quantile(val, 0.75, na.rm = TRUE),
      p95 = quantile(val, 0.95, na.rm = TRUE),
      val = ifelse(val > p95, p95, val)
    )
  
  # second trimester individual data
  t2_data <- samples %>%
    group_by(compound) %>%
    mutate(p95 = quantile(val, 0.95, na.rm = TRUE)) %>%
    filter(compound %in% comp_group & ident == id & period == "T1") %>%
    mutate(val = ifelse(val > p95, p95, val))
  
  # third trimester individual data
  t3_data <- samples %>%
    group_by(compound) %>%
    mutate(p95 = quantile(val, 0.95, na.rm = TRUE)) %>%
    filter(compound %in% comp_group & ident == id & period == "T3") %>%
    mutate(val = ifelse(val > p95, p95, val))
  
  # init plot: population distributions
  p <- ggplot(plot_data, aes(x = val)) +
    geom_histogram(alpha = 0.5)
  
  # if there is a second trimester measure, print line
  if(nrow(t2_data) > 0){
    p <- p + 
      geom_vline(
        data = t2_data, 
        mapping = aes(xintercept = val, color = "T2"), 
        lwd = 1, 
        linetype = 2
      ) 
  }
  
  # if there is a third trimester measure, print line
  if(nrow(t3_data) > 0){
    p <- p + 
      geom_vline(
        data = t3_data, 
        mapping = aes(xintercept = val, color = "T3"), 
        lwd = 1, 
        linetype = 2
      )
  }
  
  # add facet, theme and titles
  p <- p +
    facet_wrap(
      ~compound, 
      scales = "free", 
      labeller = as_labeller(comp_labels()),
      ncol = 3
    ) +
    scale_x_log10() +
    see::theme_lucid() +
    see::scale_color_material_d(palette = "complement") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      legend.title = element_blank(),
      axis.text.x = element_text(size = 8),
      axis.title = element_text(size = 10), 
      axis.text.y = element_blank(),
      legend.text = element_text(size = 10),
      strip.background = element_blank(),
    ) +
    xlab(expression(paste("niveau urinaire mesur\u00e9 en ", mu, "g/l"))) +
    ylab("")
  
  return(p)
}


