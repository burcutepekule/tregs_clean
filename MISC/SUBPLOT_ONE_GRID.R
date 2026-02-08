
variables = variables_2_plot[p_ind][[1]]

variables_name  = paste(variables, collapse = "_")

data_long = results %>%
  dplyr::select(t, tregs_on, ros_level, pat_level, rep_id, all_of(variables)) %>%
  pivot_longer(cols = all_of(variables), names_to = "variable", values_to = "value")

# At the beginning, after creating data_long
data_long = data_long %>%
  mutate(ros_level = round(ros_level, 3),
         pat_level = round(pat_level, 3))

TF_df = as.data.frame.table(TF_matrix, stringsAsFactors = FALSE)
colnames(TF_df) = c("pat_level", "ros_level", "control_percentage")

# Extract numeric values from row/column names
TF_df$pat_level = as.numeric(gsub("pat", "", TF_df$pat_level))
TF_df$ros_level = as.numeric(gsub("ros", "", TF_df$ros_level))

# Convert control_percentage to numeric (it's a proportion from 0 to 1)
TF_df$control_percentage = as.numeric(TF_df$control_percentage)

# Create a label showing the percentage
TF_df$control_label = paste0(round(TF_df$control_percentage * 100), "%")

TF_df = TF_df %>%
  mutate(ros_level = round(ros_level, 3),
         pat_level = round(pat_level, 3))

# Then filter
TF_df_filtered = TF_df %>%
  inner_join(data_long %>% distinct(ros_level, pat_level), 
             by = c("ros_level", "pat_level"))

# Create custom labeller that includes column averages
col_labeller = function(labels) {
  lapply(labels, function(x) {
    avg = col_avg$avg_pct[col_avg$ros_level == x]
    paste0("ros_level: ", x, "\n(avg: ", scales::percent(avg, accuracy = 0.1), ")")
  })
}

if(background_on[p_ind]==1){
  if(variables=='epithelial_score'){
    p = ggplot(data_long, aes(x = t, y = value)) +
      # Background with gradient coloring
      geom_rect(data = TF_df_filtered,
                aes(fill = control_percentage),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                alpha = 0.3, inherit.aes = FALSE) +
      scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                           midpoint = 0.5,
                           limits = c(0, 1),
                           name = "% Controlled",
                           labels = scales::percent) +
      # Border with gradient coloring
      geom_rect(data = TF_df_filtered,
                aes(color = control_percentage),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                fill = NA, linewidth = 1.5, inherit.aes = FALSE,
                show.legend = FALSE) +
      scale_color_gradient2(low = "darkred", mid = "orange", high = "darkgreen",
                            midpoint = 0.5,
                            limits = c(0, 1)) +
      # Add percentage text labels in top right corner of each facet
      geom_text(data = TF_df_filtered,
                aes(label = control_label),
                x = Inf, y = Inf,
                hjust = 1.1, vjust = 1.5,
                size = 4, fontface = "bold",
                color = "black",
                inherit.aes = FALSE) +
      # Reset color scale for lines
      new_scale_color() +
      # Horizontal threshold line
      geom_hline(yintercept = max_level_injury,
                 linetype = "solid",
                 color = "black",
                 linewidth = 0.5) +
      geom_hline(yintercept = epithelial_limit,
                 linetype = "solid",
                 color = "black",
                 linewidth = 0.5) +
      # Lines with agent colors
      geom_line(aes(color = variable, group = rep_id),
                alpha = alpha_plot, linewidth = 1) +
      scale_color_manual(values = agent_colors, name = "Agent") +
      # facet_grid(pat_level ~ ros_level, labeller = label_both) +
      # replace labeller to show column averages
      facet_grid(pat_level ~ ros_level, 
                 labeller = labeller(pat_level = label_both, 
                                     ros_level = col_labeller))+
      theme_minimal() +
      scale_y_log10() +
      # scale_y_log10(limits = c(NA, max_level_injury)) +
      # labs(title = title_opt, x = "Time", y = "Count")
      labs(title = "", x = "Time", y = "Count")+
      # Increase axis text size
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 18),  # top facet labels (ros_level)
        strip.text.y = element_text(size = 18),   # side facet labels (pat_level)
        legend.text = element_text(size = 14),   # legend item labels
        legend.title = element_text(size = 16)   # legend title
      )
  }else if(variables=='pathogen'){
    p = ggplot(data_long, aes(x = t, y = value)) +
      # Background with gradient coloring
      geom_rect(data = TF_df_filtered,
                aes(fill = control_percentage),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                alpha = 0.3, inherit.aes = FALSE) +
      scale_fill_gradient2(low = "red", mid = "yellow", high = "green",
                           midpoint = 0.5,
                           limits = c(0, 1),
                           name = "% Controlled",
                           labels = scales::percent) +
      # Border with gradient coloring
      geom_rect(data = TF_df_filtered,
                aes(color = control_percentage),
                xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                fill = NA, linewidth = 1.5, inherit.aes = FALSE,
                show.legend = FALSE) +
      scale_color_gradient2(low = "darkred", mid = "orange", high = "darkgreen",
                            midpoint = 0.5,
                            limits = c(0, 1)) +
      # Add percentage text labels in top right corner of each facet
      geom_text(data = TF_df_filtered,
                aes(label = control_label),
                x = Inf, y = Inf,
                hjust = 1.1, vjust = 1.5,
                size = 4, fontface = "bold",
                color = "black",
                inherit.aes = FALSE) +
      # Reset color scale for lines
      new_scale_color() +
      # Horizontal threshold line
      geom_hline(yintercept = pathogen_limit,
                 linetype = "solid",
                 color = "black",
                 linewidth = 0.5) +
      # Lines with agent colors
      geom_line(aes(color = variable, group = rep_id),
                alpha = alpha_plot, linewidth = 1) +
      scale_color_manual(values = agent_colors, name = "Agent") +
      # facet_grid(pat_level ~ ros_level, labeller = label_both) +
      # replace labeller to show column averages
      facet_grid(pat_level ~ ros_level, 
                 labeller = labeller(pat_level = label_both, 
                                     ros_level = col_labeller))+
      theme_minimal() +
      scale_y_log10() +
      # scale_y_log10(limits = c(NA, max_level_injury)) +
      # labs(title = title_opt, x = "Time", y = "Count")
      labs(title = "", x = "Time", y = "Count")+
      # Increase axis text size
      theme(
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text.x = element_text(size = 18),  # top facet labels (ros_level)
        strip.text.y = element_text(size = 18),   # side facet labels (pat_level)
        legend.text = element_text(size = 14),   # legend item labels
        legend.title = element_text(size = 16)   # legend title
      )
  }
  
}else{
  p = ggplot(data_long, aes(x = t, y = value))+
    geom_hline(yintercept = 10,
               linetype = "solid",
               color = "gray",
               linewidth = 0.5) +
    # Lines with agent colors
    geom_line(aes(color = variable, group = interaction(rep_id, variable)), 
              alpha = alpha_plot, linewidth = 1) +
    scale_color_manual(values = agent_colors, name = "Agent") +
    facet_grid(pat_level ~ ros_level, labeller = label_both) +
    theme_minimal() +
    # labs(title = title_opt, x = "Time", y = "Count")
    labs(title = "", x = "Time", y = "Count")+
    # Increase axis text size
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      strip.text.x = element_text(size = 18),  # top facet labels (ros_level)
      strip.text.y = element_text(size = 18),   # side facet labels (pat_level)
      legend.text = element_text(size = 14),   # legend item labels
      legend.title = element_text(size = 16)   # legend title
    )
}