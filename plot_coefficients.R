library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

plot_coefficients <- function(mat,
                              conf,
                              idx_range,
                              #title        = "Boxplot of Coefficients",
                              ylab         = expression(hat(beta)),
                              fill_label   = "CI Non-zero Count",
                              q_probs      = c(0.15, 0.85),
                              fill_limits  = c(0, 100),
                              fill_breaks  = seq(0, 100, by = 20),
                              viridis_opt  = "magma",
                              viridis_dir  = -1,
                              param_type   = c("beta", "omega"),
                              var_labels   = NULL) {
  param_type <- match.arg(param_type)
  
  df_long <- mat %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "Variable", values_to = "Value")
  
  df_long$Variable <- factor(df_long$Variable, levels = colnames(mat))
  
  y_stats <- quantile(df_long$Value, probs = q_probs, na.rm = TRUE)
  iqr     <- diff(y_stats)
  ymin    <- if (param_type == "omega") 0 else y_stats[1] - 1.5 * iqr
  ymax    <- y_stats[2] + 1.5 * iqr
  
  vars   <- unique(df_long$Variable)
  conf_df <- data.frame(
    Variable   = vars,
    confidence = conf[idx_range],
    stringsAsFactors = FALSE
  )
  
  df_plot <- df_long %>%
    left_join(conf_df, by = "Variable")
  
  p <- ggplot(df_plot, aes(x = Variable, y = Value, fill = confidence)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_cartesian(ylim = c(ymin, ymax)) +
    labs(
#      title = title,
      x     = "Term",
      y     = ylab,
      fill  = fill_label
    ) +
    theme(
      axis.text.x  = element_text(size=14,angle = 45, hjust = 1),
      axis.text.y  = element_text(size = 14),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 8)
    ) +
    scale_fill_viridis_c(
      option    = viridis_opt,
      direction = viridis_dir,
      limits    = fill_limits,
      breaks    = fill_breaks,
      na.value  = "grey90"
    )
  
  if (!is.null(var_labels)) {
    if (is.null(names(var_labels))) {
      names(var_labels) <- colnames(mat)
    }
    
    p <- p + scale_x_discrete(
      labels = function(x) {
        lab <- var_labels[x]
        lab[is.na(lab)] <- x[is.na(lab)]
        parse(text = lab)
      }
    )
  }
  
  p
}