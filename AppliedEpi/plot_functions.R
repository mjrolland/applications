# Plot functions
# M. Rolland
# 28/05/202

# List of compounds and their label
comp_labels <- function(){
  labs <- c("MEPA_total" = "Metylparab\u00e8ne", 
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
            "DINCHms2"   = "DINCH")
  return(labs)
}

# function to plot cohort distribution + individual levels for a given set of
# compounds
plot_hist <- function(samples, comp_group, id, my_t2, my_t3){
  
  plot_data <- samples %>%
    filter(compound %in% comp_group) %>%
    group_by(compound) %>%
    mutate(q1 = quantile(val, 0.25, na.rm = TRUE),
           q3 = quantile(val, 0.75, na.rm = TRUE),
           p95 = quantile(val, 0.95, na.rm = TRUE),
           val = ifelse(val > p95, p95, val))
  
  t2_data <- samples %>%
    group_by(compound) %>%
    mutate(p95 = quantile(val, 0.95, na.rm = TRUE)) %>%
    filter(compound %in% comp_group & ident == id & period == "T1") %>%
    mutate(val = ifelse(val > p95, p95, val))
  
  t3_data <- samples %>%
    group_by(compound) %>%
    mutate(p95 = quantile(val, 0.95, na.rm = TRUE)) %>%
    filter(compound %in% comp_group & ident == id & period == "T3") %>%
    mutate(val = ifelse(val > p95, p95, val))
  
  p <- ggplot(plot_data, aes(x = val)) +
    geom_histogram(alpha = 0.5)
  
  if(nrow(t2_data) > 0){
    p <- p +  geom_vline(data = t2_data, 
                         mapping = aes(xintercept = val, color = "T1"), 
                         lwd = 1, 
                         linetype = 2) 
  }
  
  if(nrow(t3_data) > 0){
    p <- p + geom_vline(data = t3_data, 
                        mapping = aes(xintercept = val, color = "T3"), 
                        lwd = 1, 
                        linetype = 2)
  }
  
  p <- p +
    facet_wrap(~compound, 
               scales = "free", 
               labeller = as_labeller(comp_labels()),
               ncol = 3) +
    scale_x_log10() +
    see::theme_lucid() +
    see::scale_color_material_d(palette = "complement") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      legend.title = element_blank(),
      axis.text.x = element_text(size = 8),
      #axis.title.x = element_blank(),
      axis.title = element_text(size = 10), 
      axis.text.y = element_blank(),
      legend.text = element_text(size = 10),
      strip.background = element_blank(),
    ) +
    xlab(expression(paste("niveau urinaire mesur\u00e9 en ", mu, "g/l"))) +
    ylab("")
  
  return(p)
}


plot_hist_infant <- function(samples_long, comp_group, id, my_m2, my_m12){
  
  plot_data <- samples_long %>%
    filter(compound %in% comp_group &
             !(compound == "BPS_total" & period == "M2")) %>%
    group_by(compound) %>%
    mutate(      
      p5  = quantile(val_num, 0.05, na.rm = TRUE),
      q1  = quantile(val_num, 0.25, na.rm = TRUE),
      q3  = quantile(val_num, 0.75, na.rm = TRUE),
      p95 = quantile(val_num, 0.95, na.rm = TRUE),
      val_num = case_when(
        val_num > p95 ~ p95, 
        val_num < LOD ~ LOD,
        str_detect(compound, "_ms") & val_num < p5 ~ p5, # no LOD for molar sums so we use p5 as lower limit
        TRUE ~ val_num)
    )
  
  m2_data <- plot_data %>%
    filter(ident == id & period == "M2")
  
  m12_data <- plot_data %>%
    filter(compound %in% comp_group & ident == id & period == "Y1")
  
  p <- ggplot(plot_data, 
              aes(x = val_num)) +
    geom_histogram(aes(alpha = 0.5)) +
    facet_wrap(~compound, 
               scales = "free", 
               labeller = as_labeller(comp_labels()),
               ncol = 3) +
    see::theme_lucid() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position  = "none",
      legend.title     = element_blank(),
      axis.text.x      = element_text(size = 8),
      #axis.title.x     = element_blank(),
      axis.title       = element_text(size = 10),
      # axis.text.y      = element_blank(),
      legend.text      = element_text(size = 10),
      strip.background = element_blank(),
    ) +
    scale_x_log10(name = expression(paste("niveau urinaire mesur\u00e9 en ", mu, "g/l"))) +
    scale_y_continuous(
      name   = "N enfants"
    ) +
    scale_color_manual(
      breaks = c("M2", "Y1"),
      values = c("#f44336", "#2196f3"),
    )
  
  if(nrow(m2_data) > 0){
    p <- p +  geom_vline(data     = m2_data,
                         mapping  = aes(xintercept = val_num, color = "M2"),
                         lwd      = 1,
                         linetype = 2)
  }
  
  if(nrow(m12_data) > 0){
    p <- p + geom_vline(data     = m12_data,
                        mapping  = aes(xintercept = val_num, color = "Y1"),
                        lwd      = 1,
                        linetype = 2)
  }
  
  return(p)
}

