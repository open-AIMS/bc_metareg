######################
# AUXILLIARY FUNCTIONS
######################
annotate_text <- function(label, x, y, facets = NULL, hjust = 0,
                          vjust = 0, color = "black", alpha = NA,
                          family = thm$text$family,
                          size = thm$text$size, fontface = 1, line_height = 1.0,
                          box_just = ifelse(c(x, y) < 0.5, 0, 1),
                          margin = grid::unit(size / 2, "pt"),
                          thm = theme_get()) {
# from stackoverflow
# question 22488563
# ggplot2-annotate-layer-position-in-r
  x <- scales::squish_infinite(x)
  y <- scales::squish_infinite(y)
  data <- data.frame(x = NA)
  if (!is.null(facets)) {
    data <- cbind(data, facets)
  }
  tg <- grid::textGrob(label, x = 0, y = 0, hjust = hjust, vjust = vjust,
                       gp = grid::gpar(col = alpha(color, alpha),
                                       fontsize = size, fontfamily = family,
                                       fontface = fontface,
                                       lineheight = line_height))
  ts <- grid::unit.c(grid::grobWidth(tg), grid::grobHeight(tg))
  vp <- grid::viewport(x = x, y = y, width = ts[1], height = ts[2],
                       just = box_just)
  tg <- grid::editGrob(tg, x = ts[1] * hjust, y = ts[2] * vjust, vp = vp)
  unt <- grid::unit(1, "npc") - margin * 2
  inr <- grid::grobTree(tg, vp = grid::viewport(width = unt, height = unt))
  layer(data = data, stat = StatIdentity, position = PositionIdentity,
        geom = GeomCustomAnn, inherit.aes = TRUE,
        params = list(grob = grid::grobTree(inr),
                      xmin = -Inf,
                      xmax = Inf,
                      ymin = -Inf,
                      ymax = Inf))
}

##################
# FIGURE FUNCTIONS
##################
make_fig_2 <- function(data) {
  fig_2a <- ggplot(data = data) +
    geom_density(mapping = aes(x = auto_mean, fill = Habitat), alpha = 0.5) +
    labs(x = "Autochthonous C (%)", y = "Density", fill = "Habitat") +
    scale_fill_manual(values = viridisLite::viridis(8)[c(3, 1, 5)]) +
    theme_classic()
  fig_2b <- ggplot(data = data) +
    geom_density(
      mapping = aes(x = ln_ratio_mean, fill = Habitat), alpha = 0.5
    ) +
    labs(x = "Log-ratio", y = "", fill = "Habitat") +
    scale_fill_manual(values = viridisLite::viridis(8)[c(3, 1, 5)]) +
    theme_classic()
  fig_2 <- fig_2a + fig_2b + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  ggsave("output/fig_2.pdf", fig_2, width = 6.5, height = 2.9)
  # ensure targets pipeline is triggered
  Sys.time()
}

make_fig_3 <- function(model) {
  analysis_levs <- c(
    "Bayesian mixing models", "Standard mixing models", "Other"
  )
  nd_a <- data.frame(
    Climate = c("Temperate", "Tropical"), Habitat = "Mangrove",
    Core_scaled = median(model$data$Core_scaled),
    End_m_scaled = median(model$data$End_m_scaled),
    n_scaled = median(model$data$n_scaled),
    Analysis = "Bayesian mixing models",
    ln_ratio_sd = mean(model$data$ln_ratio_sd)
  )
  fig3_a <- brms::posterior_epred(model, newdata = nd_a, re_formula = NA) |>
    data.frame() |>
    dplyr::rename(Temperate = X1, Tropical = X2) |>
    tidyr::pivot_longer(tidyselect::everything(), names_to = "Climate") |>
    ggplot(data = _) +
      geom_density(mapping = aes(x = value, fill = Climate), colour = NA,
                   alpha = 0.8, adjust = 2) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
      labs(x = "", y = "Density", fill = "", title = "by Climate") +
      scale_fill_manual(values = c("#3C5488FF", "#E64B35FF")) +
      theme_classic() +
      theme(legend.position = "top",
            legend.position.inside = c(0.8, 0.8),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 12),
            axis.title.x = element_blank())
  nd_b <- data.frame(
    Climate = "Tropical", Habitat = c("Mangrove", "Saltmarsh", "Seagrass"),
    Core_scaled = median(model$data$Core_scaled),
    End_m_scaled = median(model$data$End_m_scaled),
    n_scaled = median(model$data$n_scaled),
    Analysis = "Bayesian mixing models", ln_ratio_sd = mean(model$data$ln_ratio_sd)
  ) |>
    dplyr::mutate(
      Habitat = factor(Habitat, levels = c("Mangrove", "Saltmarsh", "Seagrass"))
    )
  fig3_b <- brms::posterior_epred(model, newdata = nd_b, re_formula = NA) |>
    data.frame() |>
    dplyr::rename(Mangrove = X1, Saltmarsh = X2, Seagrass = X3) |>
    tidyr::pivot_longer(tidyselect::everything(), names_to = "Habitat") |>
    ggplot(data = _) +
      geom_density(mapping = aes(x = value, fill = Habitat), colour = NA,
                   alpha = 0.5, adjust = 2) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
      labs(x = "", y = "", fill = "", title = "by Habitat") +
      scale_fill_manual(values = viridisLite::viridis(8)[c(3, 1, 5)]) +
      theme_classic() +
      theme(legend.position = "top",
            legend.position.inside = c(0.82, 0.8),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 12),
            axis.title.x = element_blank())
  nd_c <- data.frame(
    Climate = "Tropical", Habitat = "Mangrove",
    Core_scaled = median(model$data$Core_scaled),
    End_m_scaled = median(model$data$End_m_scaled),
    n_scaled = median(model$data$n_scaled),
    Analysis = c("Bayesian mixing models", "Standard mixing models", "Other"),
    ln_ratio_sd = mean(model$data$ln_ratio_sd)
  ) |>
    dplyr::mutate(Analysis = factor(Analysis, levels = .env$analysis_levs))
  fig3_c <- brms::posterior_epred(model, newdata = nd_c, re_formula = NA) |>
    data.frame() |>
    dplyr::rename(`Bayesian mixing models` = X1, `Standard mixing models` = X2,
                  Other = X3) |>
    tidyr::pivot_longer(tidyselect::everything(), names_to = "Analysis") |>
    ggplot(data = _) +
      geom_density(mapping = aes(x = value, colour = Analysis, fill = Analysis),
                   alpha = 0.8, adjust = 2) +
      geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
      labs(x = "", y = "", fill = "", colour = "", title = "by Method") +
      scale_colour_manual(values = c("black", "grey60", "black")) +
      scale_fill_manual(values = c("white", "grey60", "black")) +
      theme_classic() +
      theme(legend.position = "top",
            legend.position.inside = c(0.82, 0.8),
            axis.title = element_text(size = 12),
            plot.title = element_text(size = 12),
            axis.title.x = element_blank()) +
      guides(fill = guide_legend(nrow = 2, byrow = FALSE),
             colour = guide_legend(nrow = 2, byrow = FALSE))
  fig_3 <- (fig3_a + fig3_b + fig3_c) +
    patchwork::plot_annotation(tag_levels = "a",
      caption = "Ln (Auto / Allo)",
      theme = theme(
        plot.caption = element_text(size = 12, vjust = 85, hjust = 0.53)
      )
    )
  ggsave("output/fig_3.pdf", fig_3, width = 12.76, height = 3.5)
  # ensure targets pipeline is triggered
  Sys.time()
}

make_fig_4 <- function(model, mod_dat) {
  seq_vec <- function(x) {
    seq(min(x), max(x), length.out = 100)
  }
  nd_core <- data.frame(
    Climate = "Tropical", Habitat = "Mangrove",
    Analysis = "Bayesian mixing models",
    Core_scaled = seq_vec(model$data$Core_scaled),
    End_m_scaled = median(model$data$End_m_scaled),
    n_scaled = median(model$data$n_scaled),
    ln_ratio_sd = mean(model$data$ln_ratio_sd)
  )
  fig4_a <- brms::posterior_epred(model, newdata = nd_core, re_formula = NA) |>
    apply(2, ggdist::median_qi) |>
    dplyr::bind_rows() |>
    cbind(Core_scaled = seq_vec(model$data$Core_scaled)) |>
    dplyr::mutate(
      Core = Core_scaled * attr(.env$mod_dat$Core_scaled, "scaled:scale") +
        attr(.env$mod_dat$Core_scaled, "scaled:center")
    ) |>
    ggplot(data = _) +
      geom_ribbon(
        mapping = aes(x = Core, ymin = ymin, ymax = ymax), alpha = 0.5,
        fill = "grey30", colour = NA
      ) +
      geom_line(mapping = aes(x = Core, y = y)) +
      labs(x = "Core depth (m)", y = "Ln (Auto / Allo)") +
      theme_classic() +
      theme(axis.title = element_text(size = 12))
  nd_endm <- data.frame(
    Climate = "Tropical", Habitat = "Mangrove",
    Analysis = "Bayesian mixing models",
    Core_scaled = median(model$data$Core_scaled),
    End_m_scaled = seq_vec(model$data$End_m_scaled),
    n_scaled = median(model$data$n_scaled),
    ln_ratio_sd = mean(model$data$ln_ratio_sd)
  )
  fig4_b <- brms::posterior_epred(model, newdata = nd_endm, re_formula = NA) |>
    apply(2, ggdist::median_qi) |>
    dplyr::bind_rows() |>
    cbind(End_m_scaled = seq_vec(model$data$End_m_scaled)) |>
    dplyr::mutate(
      End_m = End_m_scaled * attr(.env$mod_dat$End_m_scaled, "scaled:scale") +
        attr(.env$mod_dat$End_m_scaled, "scaled:center")
    ) |>
    ggplot(data = _) +
      geom_ribbon(
        mapping = aes(x = End_m, ymin = ymin, ymax = ymax), alpha = 0.5,
        fill = "grey30", colour = NA
      ) +
      geom_line(mapping = aes(x = End_m, y = y)) +
      labs(x = "Number of end members", y = "") +
      theme_classic() +
      theme(axis.title = element_text(size = 12))
  nd_n <- data.frame(
    Climate = "Tropical", Habitat = "Mangrove",
    Analysis = "Bayesian mixing models",
    Core_scaled = median(model$data$Core_scaled),
    End_m_scaled = median(model$data$End_m_scaled),
    n_scaled = seq_vec(model$data$n_scaled),
    ln_ratio_sd = mean(model$data$ln_ratio_sd)
  )
  fig4_c <- brms::posterior_epred(model, newdata = nd_n, re_formula = NA) |>
    apply(2, ggdist::median_qi) |>
    dplyr::bind_rows() |>
    cbind(n_scaled = seq_vec(model$data$n_scaled)) |>
    dplyr::mutate(
      n = n_scaled * attr(.env$mod_dat$n_scaled, "scaled:scale") +
        attr(.env$mod_dat$n_scaled, "scaled:center")
    ) |>
    ggplot(data = _) +
      geom_ribbon(mapping = aes(x = n, ymin = ymin, ymax = ymax), alpha = 0.5,
                  fill = "grey30", colour = NA) +
      geom_line(mapping = aes(x = n, y = y)) +
      labs(x = "Number of study replicates", y = "") +
      theme_classic() +
      theme(axis.title = element_text(size = 12))
  fig_4 <- fig4_a + fig4_b + fig4_c +
    patchwork::plot_annotation(tag_levels = "a")
  ggsave("output/fig_4.pdf", fig_4, width = 7.5, height = 2.4)
  # ensure targets pipeline is triggered
  Sys.time()
}

make_fig_5 <- function(plot_data) {
  plot_data <- plot_data |>
    dplyr::filter(`Habitat Sampled` != "Unvegetated") |>
    dplyr::mutate(
      Origin = forcats::fct_recode(Origin, "Plankton" = "Phytoplankton"),
      `Habitat Sampled` = droplevels(`Habitat Sampled`)
    )
  tags <- plot_data |>
    dplyr::distinct(`Habitat Sampled`, plot_column) |>
    dplyr::arrange(plot_column) |>
    dplyr::mutate(legs = paste0("(", letters[seq_len(n())], ")"))
  labs <- tags |>
    dplyr::filter(plot_column == "a")
  my_cols <- viridisLite::viridis(8)
  names(my_cols) <- c("Saltmarsh", "Epiphytes", "Mangrove", "SPOM", "Seagrass",
                      "Plankton", "Macroalgae", "Terrestrial")
  out <- ggplot(plot_data) +
    geom_ribbon(
      mapping = aes(x = x, ymin = ymin, ymax = ymax, fill = Origin), alpha = 0.5
    ) +
    geom_line(mapping = aes(x = x, y = y, colour = Origin)) +
    labs(
      x = substitute("Contribution to soil C"["org"] * " (%)"), y = "Density",
      title = "1,000 bootstraps", fill = "Source: ", colour = "Source: "
    ) +
    scale_fill_manual(values = my_cols) +
    scale_colour_manual(values = my_cols) +
    ggh4x::facet_grid2(`Habitat Sampled` ~ plot_column, scales = "free_y",  
                       axes = "all", remove_labels = "all") +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.line = element_line(),
          strip.background = element_blank(),
          strip.text.x.top = element_blank(),
          strip.text.y.right = element_blank()) +
    guides(fill = guide_legend(nrow = 2, byrow = FALSE),
           colour = guide_legend(nrow = 2, byrow = FALSE))
  for (i in seq_len(nrow(tags))) {
    out <- out +
      annotate_text(label = unique(tags$legs)[i], facets = tags[i, ],
                    x = 0, y = 1, hjust = 0, size = 10)
  }
  for (i in seq_len(nrow(labs))) {
    out <- out +
      annotate_text(label = unique(labs$`Habitat Sampled`)[i],
                    facets = labs[i, ], x = 1, y = 1, hjust = 0, size = 10)
  }
  pys <- rep(0.4, 3); pxs <- rep(0.5, 3); x_rs <- rep(0.6, 3)
  names(pys) <- names(pxs) <- names(x_rs) <- c(
    "Mangrove", "Saltmarsh", "Seagrass"
  )
  ggsave("output/fig_5.pdf", out, width = 6.7, height = 7.5)
  # ensure targets pipeline is triggered
  Sys.time()
}

make_fig_s1 <- function(model) {
  p_pp <- brms::pp_check(model, type = "dens_overlay", ndraws = 1e3) +
    scale_colour_manual(
      values = c("black", "grey90"), labels = c("Observed", "Predicted"),
      name = ""
    ) +
    labs(x = "Study-level log-ratio mean", y = "Density") +
    theme_classic() +
    theme(legend.position = "inside",
          legend.position.inside = c(0.2, 0.8),
          legend.background = element_blank())
  p_check_data <- cbind(fitted(model, robust = FALSE), model$data) |>
    dplyr::rename(Predicted = Estimate, Observed = ln_ratio_mean)
  p_check <- ggplot(data = p_check_data) +
    geom_errorbar(
      mapping = aes(xmin = Q2.5, xmax = Q97.5, y = Observed),
      linewidth = 0.1, colour = "grey60"
    ) +
    geom_point(
      mapping = aes(x = Predicted, y = Observed, colour = Habitat,
                    shape = Analysis), size = 2
    ) +
    geom_abline(slope = 1, linetype = 2) +
    scale_colour_manual(values = viridisLite::viridis(8)[c(3, 1, 5)]) +
    labs(x = "Predicted Study-level log-ratio mean",
         y = "Observed Study-level log-ratio mean",
         colour = "Habitat:", shape = "Method:") +
    theme_classic()
  leg_colour <- ggpubr::as_ggplot(ggpubr::get_legend(p_check)[3])
  leg_shape <- ggpubr::as_ggplot(ggpubr::get_legend(p_check)[-3])
  p_check <- p_check +
    theme(legend.position = "none")
  fig_s1 <- p_pp + p_check +
    patchwork::plot_annotation(tag_levels = "a") +
    patchwork::inset_element(
      leg_colour, left = 0.02, right = 0.4, bottom = 0.6, top = 1,
      align_to = "panel", clip = FALSE, ignore_tag = TRUE
    ) +
    patchwork::inset_element(
      leg_shape, left = 0.65, right = 1, bottom = 0, top = 0.4,
      align_to = "panel", clip = FALSE, ignore_tag = TRUE
    )
  ggsave("output/fig_s1.pdf", fig_s1, width = 10.6, height = 4.5)
  # ensure targets pipeline is triggered
  Sys.time()
}

make_fig_s2 <- function(model) {
  fig_s2 <- make_brms_dharma_res(model, integerResponse = FALSE) |>
    gg_dharma(quantreg = FALSE)
  ggsave("output/fig_s2.pdf", fig_s2, width = 9, height = 4.6)
  # ensure targets pipeline is triggered
  Sys.time()
}
