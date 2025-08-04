##############
# FOR ANALYSIS
##############
empty_or_letters <- function(x) {
  l_x <- length(x)
  if (l_x == 1) {
    ""
  } else if (l_x > 1) {
    letters[seq_len(l_x)]
  }
}

truncated_sampling <- function(...) {
  x <- rnorm(1e6, ...)
  dplyr::case_when(
    x < 0 ~ 0.1, x > 100 ~ 99.9, .default = x
  ) |>
    (`/`)(100)
}

scale_w_attributes <- function(x) {
  y <- scale(x, center = TRUE, scale = TRUE)
  out <- y[, 1]
  attr(out, "scaled:center") <- attr(y, "scaled:center")
  attr(out, "scaled:scale") <- attr(y, "scaled:scale")
  out
}

filter_densities_data <- function(df_) {
  # Assume Gaussian distributions to standardise variance to SD--see report:
  #   SD = [95% CL (or CI)] * SQRT(N) / 1.96
  #   SD = [Range] / 6
  #   SD = SE * SQRT(N)
  #   SD = MAD * 1.4826
  df_ |>
    dplyr::group_by(Authors) |>
    dplyr::mutate(
      suffix = empty_or_letters(Authors), Authors = paste0(Authors, suffix)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-suffix) |>
    dplyr::select(
      `Study #`, Authors, `Habitat Sampled`, `Autochthonous`, `Autoch_var`,
      `% Mangrove`:`% Terrestrial`,
      `Variance Type`, `% MNG Var`:`% PLK Var`,
      `%saltmarsh Var`:`% Terrestrial var`, n
    ) |>
    dplyr::rename(
      `% Epiphytes` = `%Epiphytes/MPB`,
      `% Phytoplankton` = `% Plantkon`, `% Saltmarsh` = `% saltmarsh`,
      `% Mangrove Var` = `% MNG Var`,
      `% Seagrass Var` = `% SGR Var`, `% Macroalgae Var` = `% MAC Var`,
      `% Epiphytes Var` = `%EPIP Var`, `% Phytoplankton Var` = `% PLK Var`,
      `% Saltmarsh Var` = `%saltmarsh Var`, `% SPOM Var` = `%SPOM var`,
      `% Terrestrial Var` = `% Terrestrial var`
    ) |>
    tidyr::pivot_longer(
      cols = `% Mangrove`:`% Terrestrial`, names_to = "Origin",
      values_to = "Percent"
    ) |>
    tidyr::pivot_longer(cols = `% Mangrove Var`:`% Terrestrial Var`,
                        names_to = "Origin Var", values_to = "Variance") |>
    dplyr::mutate(`Origin Var` = gsub(" Var$", "", `Origin Var`)) |>
    dplyr::filter(Origin == `Origin Var`, !is.na(Percent), !is.na(Variance)) |>
    dplyr::mutate(
      SD = dplyr::case_when(
        grepl("95", `Variance Type`) ~ Variance * sqrt(n) / 1.96,
        `Variance Type` == "Range" ~ Variance / 6,
        `Variance Type` == "SE" ~ Variance * sqrt(n),
        `Variance Type` == "MAD" ~ 1.4826 * Variance,
        .default = Variance
      )
    )
}

simulate_densities <- function(niter, ...) {
  lapply(seq_len(niter), function(y, data) {
    set.seed(y)
    data |>
      dplyr::mutate(gen_per = rnorm(n(), mean = Percent, sd = SD)) |>
      split(f = ~ `Study #` + Authors + `Habitat Sampled`, drop = TRUE) |>
      lapply(function(x, y) {
        n_ <- 1
        while ((any(x$gen_per < 0) | sum(x$gen_per) > 100) & n_ <= 500) {
          set.seed(n_ * y)
          x <- dplyr::mutate(x, gen_per = rnorm(n(), mean = Percent, sd = SD))
          n_ <- n_ + 1
        }
        if (n_ > 500) {
          browser()
        }
        x
      }, y = y) |>
      dplyr::bind_rows()
  }, ...) |>
    dplyr::bind_rows(.id = "iter")
}

make_densities_plot_data <- function(sim_data) {
  col_1 <- c("Macroalgae", "Mangrove", "Seagrass", "Saltmarsh")
  plot_data <- sim_data |>
    split(f = ~ `Habitat Sampled` + Origin + iter, drop = TRUE) |>
    purrr::map(function(x) {
      out <- density(x$gen_per, n = 512, from = 0, to = 100, adjust = 2)
      if (!any(out$y > 0.1)) {
        data.frame(
          `Habitat Sampled` = x$`Habitat Sampled`[1], Origin = x$Origin[1],
          iter = x$iter[1], x = out$x, y = out$y, check.names = FALSE
        )
      }
    }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      Origin = gsub("% ", "", Origin, fixed = TRUE), iter = as.factor(iter),
      plot_column = ifelse(Origin %in% .env$col_1, "a", "b"),
      `Habitat Sampled` = forcats::fct_relevel(`Habitat Sampled`,
        c("Saltmarsh", "Mangrove", "Seagrass", "Unvegetated")
      ),
      Origin = forcats::fct_relevel(Origin,
        c("Saltmarsh", "Mangrove", "Seagrass", "Macroalgae", "Epiphytes",
          "SPOM", "Phytoplankton", "Terrestrial")
      )
    )
  plot_data |>
    dplyr::group_by(`Habitat Sampled`, Origin, plot_column, x) |>
    dplyr::summarise(ggdist::mean_hdci(y)) |>
    dplyr::ungroup()
}

###########
# MODELLING
###########
filter_model_data <- function(dat_raw) {
  analysis_levs <- c(
    "Bayesian mixing models", "Standard mixing models", "Other"
  )
  dat <- dat_raw |>
    dplyr::filter(is.na(`Excluded?`)) |>
    dplyr::select(
      `Study #`, Authors, Climate, `Habitat Sampled`, n,
      Core = `Core depth midpoint (m)`, End_m = `End-members attempted`,
      Analysis = `Analysis method`, Location, `% Mangrove`:`% Terrestrial`,
      `Variance Type`, `% MNG Var`:`% PLK Var`,
      `%saltmarsh Var`:`% Terrestrial var`, Autochthonous, Allochthonous
    ) |>
    dplyr::rename(Study = `Study #`, Habitat = `Habitat Sampled`) |>
    dplyr::filter(
      Habitat != "Unvegetated", !is.na(Core),
      !is.na(End_m), !is.na(Analysis), !is.na(Climate)
    ) |>
    dplyr::group_by(Authors) |>
    dplyr::mutate(
      suffix = empty_or_letters(Authors), Authors = paste0(Authors, suffix)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-suffix)
  var_adj <- dat |>
    dplyr::select(
      Study, Authors, Habitat, `% Mangrove`:`% Terrestrial`, `Variance Type`,
      `% MNG Var`:`% PLK Var`, `%saltmarsh Var`:`% Terrestrial var`, n
    ) |>
    dplyr::group_by(Study, Authors, Habitat) |>
    dplyr::ungroup() |>
    dplyr::rename(
      `% Epiphytes` = `%Epiphytes/MPB`,
      `% Saltmarsh` = `% saltmarsh`,
      `% Mangrove Var` = `% MNG Var`,
      `% Seagrass Var` = `% SGR Var`, `% Macroalgae Var` = `% MAC Var`,
      `% Epiphytes Var` = `%EPIP Var`, `% Phytoplankton Var` = `% PLK Var`,
      `% Saltmarsh Var` = `%saltmarsh Var`, `% SPOM Var` = `%SPOM var`,
      `% Terrestrial Var` = `% Terrestrial var`
    ) |>
    tidyr::pivot_longer(
      cols = `% Mangrove`:`% Terrestrial`, names_to = "Origin",
      values_to = "Percent"
    ) |>
    tidyr::pivot_longer(cols = `% Mangrove Var`:`% Terrestrial Var`,
                        names_to = "Origin Var", values_to = "Variance") |>
    dplyr::mutate(`Origin Var` = gsub(" Var$", "", `Origin Var`)) |>
    dplyr::filter(Origin == `Origin Var`, !is.na(Percent), !is.na(Variance)) |>
    dplyr::mutate(
      SD = dplyr::case_when(
        grepl("95", `Variance Type`) ~ Variance * sqrt(n) / 1.96,
        `Variance Type` == "Range" ~ Variance / 6,
        `Variance Type` == "SE" ~ Variance * sqrt(n),
        `Variance Type` == "MAD" ~ 1.4826 * Variance,
        .default = Variance
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(`Origin Var` = str_trim(str_remove(`Origin Var`, "%"))) |>
    dplyr::mutate(
      Source = dplyr::case_when(
        str_detect(Habitat, `Origin Var`) == TRUE ~ "auto", .default = "allo"
      )
    ) |>
    dplyr::group_by(Study, Authors, Habitat, Source) |>
    dplyr::summarise(mean = sum(Percent), sd = sqrt(sum(SD ^ 2))) |>
    tidyr::pivot_wider(
      names_from = Source, names_glue = "{Source}_{.value}",
      values_from = c(mean, sd)
    ) |>
    dplyr::mutate(allo_mean = 100 - auto_mean) |>
    dplyr::select(-allo_sd) |>
    tidyr::drop_na(auto_mean)
  dat |>
    dplyr::left_join(var_adj) |>
    dplyr::mutate(
      Analysis = dplyr::case_match(Analysis,
        "Mass balance model" ~ "Standard mixing models",
        "Linear mixing model" ~ "Standard mixing models",
        "Bayesian Mixing Model" ~ "Bayesian mixing models",
        "Bayesian Mixing Model (SimmR)" ~ "Bayesian mixing models",
        "Biomarker rel. abund./ratio" ~ "Other", "Mixing diagram" ~ "Other",
        "Not conducted/described" ~ "Other", .default = Analysis
      )
    ) |>
    tidyr::drop_na(auto_mean) |>
    dplyr::group_by(Study, Authors, Habitat) |>
    tidyr::nest() |>
    dplyr::mutate(
      nested = purrr::map(data,
        function(x) {
          out_v <- purrr::pmap(list(x$auto_mean, x$auto_sd), function(v1, v2) {
            truncated_sampling(mean = v1, sd = v2)
          }) |>
            unlist()
          out_z <- 1 - out_v
          if (any(is.na(log(out_v / out_z)))) browser() # is there an issue?
          out <- log(out_v / out_z) # log-ratio allo / auto
          list(ln_ratio_mean = mean(out), ln_ratio_sd = sd(out))
        }
      )
    ) |>
    tidyr::unnest_wider(nested) |>
    tidyr::unnest_wider(data) |>
    dplyr::ungroup() |>
    dplyr::select(Study:Location, allo_mean:ln_ratio_sd) |>
    dplyr::mutate(
      ln_ratio_sd = dplyr::case_match(ln_ratio_sd,
        0 ~ 0.0001, .default = ln_ratio_sd
      )
    ) |>
    dplyr::mutate(
      Habitat = factor(
        Habitat, levels = c("Mangrove", "Saltmarsh", "Seagrass")
      ),
      Analysis = factor(Analysis, levels = .env$analysis_levs),
      Climate = factor(Climate, levels = c("Tropical", "Temperate"))
    ) |>
    tidyr::drop_na() |>
    dplyr::mutate(
      dplyr::across(c(n, Core, End_m), identity, .names = "{.col}_scaled"),
      dplyr::across(ends_with("_scaled"), scale_w_attributes)
    )
}

fit_metreg <- function(data) {
  priors <- brms::prior(normal(0, 1), class = "b") +
    brms::prior(normal(0.5, 0.1), class = "sd") +
    brms::prior(normal(0.5, 0.1), class = "sigma")
  model <- brms::brm(
    brms::bf(ln_ratio_mean | se(ln_ratio_sd, sigma = TRUE) ~ Climate + Habitat +
               Core_scaled + End_m_scaled + n_scaled + Analysis + (1 | Authors),
             center = FALSE), data = data, family = gaussian(),
    prior = priors, iter = 1e4, warmup = 5e3, cores = 4, chains = 4,
    sample_prior = "yes", seed = 10, save_pars = save_pars(all = TRUE),
    control = list(max_treedepth = 20, adapt_delta = 0.9)
  )
  brms::add_criterion(model, "loo")
}

fit_null_metreg <- function(data) {
  # Null model
  priors_null <- brms::prior(normal(0, 1), class = "b") +
    brms::prior(normal(0.5, 0.1), class = "sd") +
    brms::prior(normal(0.5, 0.1), class = "sigma")
  b0 <- brms::brm(
    brms::bf(ln_ratio_mean | se(ln_ratio_sd, sigma = TRUE) ~ 1 + (1 | Authors),
             center = FALSE),
    data = data, family = gaussian(), prior = priors_null, iter = 1e4,
    warmup = 5e3, cores = 4, chains = 4, sample_prior = "yes", seed = 10,
    save_pars = save_pars(all = TRUE),
    control = list(max_treedepth = 20, adapt_delta = 0.9)
  )
  brms::add_criterion(b0, "loo")
}

###############################
###############################
#### ggplot2-style DHARMa plots
###############################
###############################
make_brms_dharma_res <- function(brms_model, seed = 10, ...) {
  # equivalent to `simulateResiduals(lme4_model, use.u = FALSE)`
  # cores are set to 1 just to ensure reproducibility
  options(mc.cores = 1)
  on.exit(options(mc.cores = parallel::detectCores()))
  response <- brms::standata(brms_model)$Y
  ndraws <- nrow(as_draws_df(brms_model))
  manual_preds_brms <- matrix(0, ndraws, nrow(brms_model$data))
  random_terms <- insight::find_random(
    brms_model, split_nested = TRUE, flatten = TRUE
  )
  # for this to have a similar output to `glmmTMB`'s default, we need to
  #   create new levels in the hierarchical variables, so then we can
  #   use `allow_new_levels = TRUE` and `sample_new_levels = "gaussian"` in
  #   `brms::posterior_epred`. This is equivalent to
  #   `simulateResiduals(lme4_model, use.u = FALSE)`. See details in
  #   `lme4:::simulate.merMod` and `glmmTMB:::simulate.glmmTMB`
  new_data <- brms_model$data |>
    dplyr::mutate(across(
      all_of(random_terms), \(x)paste0("NEW_", x) |> as.factor()
    ))
  set.seed(seed)
  brms_sims <- brms::posterior_predict(
    brms_model, re_formula = NULL, newdata = new_data,
    allow_new_levels = TRUE, sample_new_levels = "gaussian"
  ) |>
    t()
  fitted_median_brms <- apply(brms_sims, 1, median)
  DHARMa::createDHARMa(
    simulatedResponse = brms_sims,
    observedResponse = response,
    fittedPredictedResponse = fitted_median_brms,
    ...
  )
}

check_dots <- DHARMa:::checkDots
dispersion_fct <- DHARMa:::testDispersion
ensure_dharma <- DHARMa:::ensureDHARMa
ensure_predictor <- DHARMa:::ensurePredictor
get_p_val <- DHARMa:::getP
levene_test <- DHARMa:::leveneTest_formula
outliers_fct <- DHARMa:::testOutliers
uniformity_fct <- DHARMa:::testUniformity
test_categorical <- DHARMa:::testCategorical
test_quantiles <- DHARMa:::testQuantiles

plot_qq_unif <- function(sim_listrix, test_uniformity = TRUE,
                         test_outliers = TRUE, test_dispersion = TRUE, ...) {
  sim_listrix <- ensure_dharma(sim_listrix, convert = "Model")
  res_ <- data.frame(y = sim_listrix$scaledResiduals) %>%
    dplyr::mutate(x = seq_len(dplyr::n()) / (dplyr::n() + 1))
  res_ <- qqplot(res_$x, res_$y, plot.it = FALSE) %>%
    as.data.frame
  p_ <- ggplot(data = res_) +
    geom_point(mapping = aes(x = x, y = y), shape = 16,
               colour = "skyblue", alpha = 0.8) +
    geom_abline(slope = 1, linetype = 2) +
    labs(x = "Expected", y = "Observed", subtitle = "QQ plot residuals") +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    theme_bw() +
    theme(plot.title = element_text(size = 9))
  if (test_uniformity) {
    tmp_ks <- uniformity_fct(sim_listrix, plot = FALSE)
    k_lab_a <- paste("KS test: p =", round(tmp_ks$p.value, digits = 5))
    k_lab_b <- paste("Deviation:", ifelse(tmp_ks$p.value < 0.05, "Significant",
                                          "N.S."))
    p_ <- p_ +
      annotate("text", x = 0, y = 0.98, hjust = 0, vjust = 0.5, label = k_lab_a,
               colour = ifelse(tmp_ks$p.value < 0.05, "tomato", "black"),
               size = 3) +
      annotate("text", x = 0, y = 0.90, hjust = 0, vjust = 0.5, label = k_lab_b,
               colour = ifelse(tmp_ks$p.value < 0.05, "tomato", "black"),
               size = 3)
  }
  if (test_outliers) {
    tmp_ot <- outliers_fct(sim_listrix, plot = FALSE)
    ot_lab_a <- paste("Outlier test: p =", round(tmp_ot$p.value, digits = 5))
    ot_lab_b <- paste("Deviation:", ifelse(tmp_ot$p.value < 0.05, "Significant",
                                           "N.S."))
    p_ <- p_ +
      annotate("text", x = 0, y = 0.78, hjust = 0, vjust = 0.5, size = 3,
               label = ot_lab_a, colour = ifelse(tmp_ot$p.value < 0.05,
                                                 "tomato", "black")) +
      annotate("text", x = 0, y = 0.70, hjust = 0, vjust = 0.5, size = 3,
               label = ot_lab_b, colour = ifelse(tmp_ot$p.value < 0.05,
                                                 "tomato", "black"))
    }
  if (test_dispersion) {
    tmp_ds <- dispersion_fct(sim_listrix, alternative = "two.sided",
                             plot = FALSE)
    ds_lab_a <- paste("Dispersion test: p =", round(tmp_ds$p.value, digits = 5))
    ds_lab_b <- paste("Deviation:", ifelse(tmp_ds$p.value < 0.05, "significant",
                                           "N.S."))
    p_ <- p_ +
      annotate("text", x = 1, y = 0.12, hjust = 1, vjust = 0.5, size = 3,
               label = ds_lab_a, colour = ifelse(tmp_ds$p.value < 0.05,
                                                 "tomato", "black")) +
      annotate("text", x = 1, y = 0.04, hjust = 1, vjust = 0.5, size = 3,
               label = ds_lab_b, colour = ifelse(tmp_ds$p.value < 0.05,
                                                 "tomato", "black"))
  }
  p_
}

plot_residuals <- function(sim_list, form = NULL, quantreg = NULL, rank = TRUE,
                           as_factor = NULL, smooth_scatter = NULL,
                           quantiles = c(0.25, 0.5, 0.75), ...) {
  ##### Checks #####
  if (is.null(form)) {
    xlab_ <- check_dots("xlab", ifelse(rank, "Model predictions (rank transformed)", "Model predictions"), ...)
  } else {
    xlab_ <- "Predictor"
  }
  sim_list <- ensure_dharma(sim_list, convert = TRUE)
  res <- sim_list$scaledResiduals
  if (inherits(form, "DHARMa")) {
    stop("DHARMa::plot_residuals > argument form cannot be of class DHARMa. ",
         "Note that the syntax of plot_residuals has changed since DHARMa ",
         "0.3.0. See ?plot_residuals.")
  }
  pred <- ensure_predictor(sim_list, form)

  ##### Rank transform and factor conversion #####
  if (!is.factor(pred)) {
    if (rank) {
      pred <- rank(pred, ties.method = "average")
      pred <- pred / max(pred)
    }
    nuniq <- length(unique(pred))
    ndata <- length(pred)
    if (is.null(as_factor)) {
      as_factor <- (nuniq == 1) | (nuniq < 10 & ndata / nuniq > 10)
    }
    if (as_factor) {
      pred <- factor(pred)
    }
  }

  ##### Residual scatter plots #####
  if (is.null(quantreg)) {
    if (length(res) > 2000) {
      quantreg <- FALSE
    } else {
       quantreg <- TRUE
    }
  }
  switch_scatter <- 1e4
  if (is.null(smooth_scatter)) {
    if (length(res) > switch_scatter) {
      smooth_scatter <- TRUE
    } else {
      smooth_scatter <- FALSE
    }
  }
  blackcol <- rgb(0, 0, 0,
                  alpha = max(0.1, 1 - 3 * length(res) / switch_scatter))
  if (is.factor(pred)) {
    out <- cat_gg_res(sim_list = sim_list, cat_pred = pred,
                      quantiles = quantiles)
  } else if (smooth_scatter) {
    def_col_ <- ifelse(res == 0 | res == 1, 2, blackcol)
    df_ <- data.frame(x = pred, y = res, col = def_col_ == 2)
    out <- ggplot(data = df_, mapping = aes(x = x, y = y)) +
      stat_density2d(mapping = aes(fill = ..density..^0.25),
                     geom = "tile", contour = FALSE, n = 200) +
      geom_point(mapping = aes(x = x, y = y), size = 0.2, shape = 16,
                 colour = "black") +
      geom_point(data = df_ %>% dplyr::filter(col),
                 mapping = aes(x = x, y = y), size = 0.2, shape = 16,
                 colour = "tomato") +
      scale_fill_continuous(low = "white", high = "dodgerblue4") +
      ylim(c(0, 1)) +
      theme_bw()
  } else {
    symbol <- ifelse(res == 0 | res == 1, "Yes", "No")
    df_ <- data.frame(x = pred, y = res, extreme = symbol)
    out <- ggplot() +
      geom_point(data = df_ %>% dplyr::filter(extreme == "No"),
                 mapping = aes(x = x, y = y), colour = "black", shape = 1,
                 size = 1) +
      geom_point(data = df_ %>% dplyr::filter(extreme == "Yes"),
                 mapping = aes(x = x, y = y), colour = "tomato", shape = 8,
                 size = 1) +
      ylim(c(0, 1)) +
      labs(y = "Standardized residual", x = xlab_) +
      theme_bw()
  }

  # ##### Quantile regressions #####
  tit_ <- check_dots("main", "Residual vs. predicted", ...)
  tmp_ <- NULL
  if (is.numeric(pred)) {
    if (!quantreg) {
      out <- out +
        labs(title = tit_) +
        geom_hline(yintercept = quantiles, linetype = 2, colour = "lightgrey")
      sspl_ <- try({
        smooth.spline(pred, res, df = 10)
      }, silent = TRUE)
      if (inherits(sspl_, "try-error")) {
        out <- out +
          geom_line(data = data.frame(x = sspl_$x, y = sspl_$y),
                    mapping = aes(x = x, y = y), linetype = 2,
                    colour = "tomato") +
          geom_hline(yintercept = 0.5, linetype = 1, colour = "tomato")
      }
    } else {
      tmp_ <- test_quantiles(sim_list, pred, quantiles = quantiles,
                             plot = FALSE)
      if (any(tmp_$pvals < 0.05, na.rm = TRUE)) {
        tit_ <- paste(tit_, "Quantile deviations detected (red curves)",
                      sep = "\n")
        if (tmp_$p.value <= 0.05) {
          tit_ <- paste(tit_, "Combined adjusted quantile test: Significant",
                        sep = "\n")
        } else {
          tit_ <- paste(tit_, "Combined adjusted quantile test: N.S.",
                        sep = "\n")
        }
        maincol <- "tomato"
      } else {
        tit_ <- paste(tit_, "No significant problems detected", sep = "\n")
        maincol <- "black"
      }
      out <- out +
        labs(title = tit_) +
        theme(plot.title = element_text(size = 9, colour = maincol))
      for (i in seq_along(quantiles)) {
        line_col <- ifelse(tmp_$pvals[i] <= 0.05 & !(is.na(tmp_$pvals[i])),
                           "tomato", "black")
        pol_df_ <- data.frame(
          x = c(tmp_$predictions$pred, rev(tmp_$predictions$pred)),
          y = c(tmp_$predictions[, 2 * i] - tmp_$predictions[, 2 * i + 1],
                rev(tmp_$predictions[, 2 * i] + tmp_$predictions[, 2 * i + 1]))
        )
        pred_df_ <- data.frame(x = tmp_$predictions$pred,
                               y = tmp_$predictions[, 2 * i])
        out <- out +
          geom_hline(yintercept = quantiles[i], colour = line_col, lwd = 0.5,
                     linetype = 2) +
          geom_polygon(data = pol_df_, mapping = aes(x = x, y = y),
                       fill = "#00000020") +
          geom_line(data = pred_df_, mapping = aes(x = x, y = y),
                    colour = line_col, lwd = 0.5)        
      }
      subt_ <- c(paste("Quantile test: p =", round(tmp_$p.value, digits = 5)),
                 paste("Deviation:", ifelse(tmp_$p.value < 0.05, "Significant",
                                             "N.S.")))
      subt_col_ <- ifelse(tmp_$p.value < 0.05, "tomato", "black")
      out <- out +
        labs(subtitle = subt_) +
        theme(plot.subtitle = element_text(size = 9, colour = subt_col_))
    }
  }
  out
}

cat_gg_res <- function(sim_list, cat_pred, quantiles = c(0.25, 0.5, 0.75),
                       plot = TRUE) {
  sim_list <- ensure_dharma(sim_list, convert = TRUE)
  cat_pred <- as.factor(cat_pred)
  out <- list()
  out$uniformity$details <- suppressWarnings(by(sim_list$scaledResiduals,
                                                cat_pred, ks.test, "punif",
                                                simplify = TRUE))
  out$uniformity$p.value <- rep(NA, nlevels(cat_pred))
  for (i in seq_len(nlevels(cat_pred))) {
    out$uniformity$p.value[i] <- out$uniformity$details[[i]]$p.value
  }
  out$uniformity$p.value.cor <- p.adjust(out$uniformity$p.value)
  if (nlevels(cat_pred) > 1) {
    out$homogeneity <- levene_test(sim_list$scaledResiduals ~ cat_pred)
  }
  if (plot) {
    if (length(out) > 1) {
        subt_ <- ifelse(out$homogeneity$`Pr(>F)`[1] < 0.05,
                        "Levene Test for homongeneity of variance: significant",
                        "Levene Test for homongeneity of variance: N.S.")
        subt_col_ <- ifelse(out$homogeneity$`Pr(>F)`[1] < 0.05, "red", "black")
    } else {
      subt_ <- ""
      subt_col_ <- "black"
    }
    tit_ <- ifelse(any(out$uniformity$p.value.cor < 0.05),
                   "Within-group deviations from uniformity: significant (red)",
                   "Within-group deviation from uniformity: N.S.")
    tit_col_ <- ifelse(any(out$uniformity$p.value.cor < 0.05), "red", "black")
    p_ <- data.frame(y = sim_list$scaledResiduals, x = cat_pred) %>%
      ggplot(data = .) +
        geom_boxplot(mapping = aes(x = x, y = y),
                     colour = ifelse(out$uniformity$p.value.cor < 0.05, "red",
                                     "black")) +
        theme_bw() +
        labs(x = "Categorical predictor", y = "Scaled residuals",
             title = tit_, subtitle = subt_) +
        scale_y_continuous(breaks = c(0, quantiles, 1), limits = c(0, 1)) +
        geom_hline(yintercept = quantiles, linetype = 2, colour = "lightgrey") +
        theme(plot.title = element_text(size = 9, colour = tit_col_),
              plot.subtitle = element_text(size = 9, colour = subt_col_))
  }
  p_
}

gg_disp_hist <- function(sim_list, alternative = c("two.sided", "greater",
                                                   "less")) {
  sim_list <- ensure_dharma(sim_list, convert = "Model")
  expected_var <- sd(sim_list$simulatedResponse) ^ 2
  spread <- function(x, sim_list, expected_var) {
    var(x - sim_list$fittedPredictedResponse) / expected_var
  }
  alternative <- match.arg(alternative)
  observed <- spread(sim_list$observedResponse, sim_list, expected_var)
  simulated <- apply(sim_list$simulatedResponse, 2, spread, sim_list,
                     expected_var)
  p_val_ <- get_p_val(simulated = simulated, observed = observed,
                      alternative = alternative)
  x_lab_ <- paste("Simulated values, red line = fitted model. p-value (",
                  alternative, ") = ", p_val_, sep = "")
  data.frame(simulated = simulated) %>%
    ggplot(data = .) +
      geom_histogram(mapping = aes(x = simulated), colour = "grey60",
                     bins = max(round(sim_list$nSim / 5), 20)) +
      geom_vline(xintercept = observed, colour = "tomato") +
      labs(x = x_lab_, y = "",
           title = "DHARMa nonparametric dispersion test via sd of residuals",
           subtitle = "Fitted vs. Simulated") +
      theme_bw() +
      theme(plot.title = element_text(size = 9),
            plot.subtitle = element_text(size = 9))
}

# Brendans edited dispersion histogram function (included cutoff_hdci so that we can still see the distribution of simulated values when the spread is very large but the vast majority of values is continaed in a small range)
gg_dispersion_hist <- function(
  sim_list, 
  alternative = c("two.sided", "greater", "less"), 
  cutoff_hdci = NULL,  # Only plot simulated values within the specified highest continuous density interval (HDCI); decimal between 0 and 1, e.g. cutoff_hdci = 0.99 plots simulated values within the 99% HDCI. 
  wrap_subtitle = NULL # Integer; if specified, lines will be wrapped to the given number of characters
) {
  sim_list <- ensure_dharma(sim_list, convert = "Model")
  expected_var <- sd(sim_list$simulatedResponse) ^ 2
  spread <- function(x, sim_list, expected_var) {
    var(x - sim_list$fittedPredictedResponse) / expected_var
  }
  alternative <- match.arg(alternative)
  observed <- spread(sim_list$observedResponse, sim_list, expected_var)
  simulated <- apply(sim_list$simulatedResponse, 2, spread, sim_list, expected_var)
  p_val_ <- get_p_val(simulated = simulated, observed = observed, alternative = alternative)
  x_lab_ <- paste0("Simulated values, red line = fitted model. p-value (", alternative, ") = ", p_val_)
  subtitle <- "Fitted vs. Simulated"
  if (!is.null(cutoff_hdci)) {
    cutoff_vals <- ggdist::hdci(simulated, .width = cutoff_hdci)
    simulated <- simulated[simulated > cutoff_vals[1] & 
      simulated < cutoff_vals[2]]
    subtitle <- paste0(subtitle, ". Only showing simulated values within the ", cutoff_hdci*100, "% highest continuous density interval (HDCI).")
  }
  if(!is.null(wrap_subtitle)) subtitle <- stringr::str_wrap(subtitle, wrap_subtitle)
  nbins <- max(floor(length(simulated)/5), 20)  
  simulated |> 
    data.frame() |> 
    ggplot(aes(simulated)) +
    geom_histogram(color = "black", fill = "black", bins = nbins) +
    geom_vline(
      xintercept = observed, 
      colour = "red",
      linetype = "longdash", 
      linewidth = 0.8
    ) +
    labs(
      x = x_lab_, 
      y = "", 
      title = "DHARMa nonparametric dispersion test via sd of residuals", 
      subtitle = subtitle
    )
}

gg_zero_inflation_hist <- function(sim_list, alternative = c("two.sided", "greater", "less"), wrap_subtitle = NULL) {
# NOTES: Zero inflation is a property of the model, not the data. It occurs when the model does not appropriately capture the data-generating process leading to the number of zeros in the observed data. This test considers the distribution of the number of zeros in many simulations of data from the fitted model. If the number of zeros in the observed data does not conform to the distribution of zeros from simulated datasets, we have evidence of zero-inflation (i.e. the fitted model is misspecified). This test is likely better suited for detecting zero-inflation than the standard DHARMa residual/QQ plots, but note that overdispersion can also lead to excess zeros. Therefore, seeing too many zeros is not a reliable diagnostics for moving towards a zero-inflated model. A reliable differentiation between overdispersion and zero-inflation will usually only be possible when directly comparing alternative models (e.g. through residual comparison/model selection of a model with and without zero-inflation, or by simply fitting a model with zero-inflation and looking at the parameter estimate for the zero-inflation)
# ARGUMENTS: 
#  - sim_list: an object made with either DHARMa::simulateResiduals() or make_brms_dharma_res; the simulated residuals
#  - alternative: string; the type of significance test to run (defaults to two.sided)
#  - wrap_subtitle: integer; if specified, lines will be wrapped to the given number of characters
  require(ggplot2)
  require(stringr)
  require(DHARMa)
  if (length(alternative) == 3) {
    message("`alternative` argument not provided, using `alternative = 'two.sided'`")
    alternative <- "two.sided"
  }
  count_zeros <- function(x) sum(x == 0)
  sim_list <- ensure_dharma(sim_list, convert = "Model")
  observed <- count_zeros(sim_list$observedResponse)
  simulated <- apply(sim_list$simulatedResponse, 2, count_zeros)
  p <- get_p_val(simulated = simulated, observed = observed, alternative = alternative)
  alt_formatted <- alternative |> 
    stringr::str_replace("\\.", "-") |> 
    stringr::str_to_sentence()
  alt_formatted <- ifelse(alternative == "two.sided", alt_formatted, paste0("One-sided (", alternative, ")"))
  subtitle <- paste0(
    "Number of zeros in observed data versus in ", sim_list$nSim, 
    " simulations of equivalent sample size (n = ", sim_list$nObs, 
    ") from the fitted model.\n", alt_formatted, " p-value = ", signif(p, 3), 
    " for H0: fitted model is a reasonable representation of the Y = 0 data generating process."
  )
  if(!is.null(wrap_subtitle)) subtitle <- stringr::str_wrap(subtitle, wrap_subtitle)
  data.frame(simulated) |> 
    ggplot(aes(simulated)) + 
    geom_histogram(
      stat = "count",
      binwidth = 1,
      fill = "black"
    ) + 
    geom_vline(
      xintercept = observed, 
      color = "red", 
      linetype = "longdash", 
      linewidth = 0.8
    ) + 
    labs(
      title = "DHARMa non-parametric zero-inflation test",
      subtitle = subtitle,
      y = "Frequency", 
      x = "Number of zeros in simulated data\n(observed number of zeros shown as red dashed line)"
    )
}

gg_dharma <- function(x, title = "DHARMa residual diagnostics", ...) {
  require(patchwork)
  a_ <- plot_qq_unif(x)
  b_ <- plot_residuals(x, ...)
  a_ + b_ + patchwork::plot_annotation(title = title)
}
