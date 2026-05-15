utils::globalVariables(c("xmid", "ymid", "sz", "x", "ylo", "yhi",
                          "est", "label", "lo", "hi"))

plot.rdrobust = function(x, y, x_run, nbins = 20,
                         ci          = TRUE,
                         show_effect = FALSE,
                         title       = NULL,
                         x.label     = "Running Variable",
                         y.label     = "Outcome",
                         x.lim       = NULL,
                         y.lim       = NULL,
                         col.l       = "#3B7DD8",
                         col.r       = "#D95F3B",
                         base_size   = 14,
                         ...) {

  obj <- x   # rdrobust object; rename for clarity

  # --- extract key quantities --------------------------------------------------
  c_val   <- obj$c
  h_l     <- obj$bws["h", "left"]
  h_r     <- obj$bws["h", "right"]
  coefs_l <- as.numeric(obj$beta_Y_p_l)
  coefs_r <- as.numeric(obj$beta_Y_p_r)
  V_l     <- obj$V_cl_l
  V_r     <- obj$V_cl_r
  kernel  <- obj$kernel

  tau_cl   <- obj$coef["Conventional",  1]
  ci_rb    <- obj$ci["Robust", ]
  lev      <- obj$level
  z_crit   <- qnorm(1 - (1 - lev / 100) / 2)
  is_fuzzy <- !is.null(obj$tau_T)

  # ggplot text sizes: annotate/geom_text use mm, ~1pt = 0.353mm
  txt_main <- base_size * 0.353        # body annotations
  txt_sub  <- (base_size - 2) * 0.353  # footnote / secondary

  # --- Healy-inspired theme ----------------------------------------------------
  theme_rd <- function() {
    ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        panel.border        = ggplot2::element_blank(),
        axis.line           = ggplot2::element_line(color = "grey30",
                                                    linewidth = 0.4),
        panel.grid.major    = ggplot2::element_line(color = "grey92",
                                                    linewidth = 0.4),
        panel.grid.minor    = ggplot2::element_blank(),
        plot.title          = ggplot2::element_text(
                                face = "bold",
                                size = base_size + 2,
                                margin = ggplot2::margin(b = 6)),
        axis.title          = ggplot2::element_text(size = base_size),
        axis.text           = ggplot2::element_text(size = base_size - 1,
                                                    color = "grey30"),
        plot.margin         = ggplot2::margin(12, 16, 8, 8),
        plot.caption        = ggplot2::element_text(size = base_size - 3,
                                                    color = "grey50",
                                                    hjust = 0,
                                                    margin = ggplot2::margin(t = 6))
      )
  }

  # --- kernel weight function --------------------------------------------------
  kern_w <- function(u) {
    u <- abs(u)
    switch(tolower(substr(kernel, 1, 3)),
      "tri" = pmax(1 - u,            0),
      "epa" = pmax(0.75 * (1 - u^2), 0),
      "uni" = as.numeric(u <= 1) * 0.5,
      pmax(1 - u, 0)
    )
  }

  # --- filter to bandwidth window + drop NAs -----------------------------------
  keep   <- !is.na(x_run) & !is.na(y) &
            (x_run >= c_val - h_l) & (x_run <= c_val + h_r)
  x_w    <- x_run[keep]
  y_w    <- y[keep]
  side   <- ifelse(x_w < c_val, "left", "right")

  kw <- ifelse(side == "left",
               kern_w((x_w - c_val) / h_l),
               kern_w((x_w - c_val) / h_r))

  # --- binned means + average kernel weight per bin ----------------------------
  make_bins <- function(xv, yv, wv, nb) {
    brks <- seq(min(xv), max(xv), length.out = nb + 1)
    grp  <- cut(xv, brks, include.lowest = TRUE)
    out  <- data.frame(xmid = as.numeric(tapply(xv, grp, mean, na.rm = TRUE)),
                       ymid = as.numeric(tapply(yv, grp, mean, na.rm = TRUE)),
                       wavg = as.numeric(tapply(wv, grp, mean, na.rm = TRUE)))
    out[stats::complete.cases(out), ]
  }

  idx_l  <- side == "left"
  idx_r  <- side == "right"
  bins_l <- make_bins(x_w[idx_l], y_w[idx_l], kw[idx_l], nbins)
  bins_r <- make_bins(x_w[idx_r], y_w[idx_r], kw[idx_r], nbins)

  # rescale kernel weights to point size range [1.5, 5]
  all_w  <- c(bins_l$wavg, bins_r$wavg)
  w_min  <- min(all_w, na.rm = TRUE)
  w_rng  <- max(all_w, na.rm = TRUE) - w_min
  if (w_rng == 0) w_rng <- 1
  rescale_sz <- function(w) 1.5 + 3.5 * (w - w_min) / w_rng
  bins_l$sz  <- rescale_sz(bins_l$wavg)
  bins_r$sz  <- rescale_sz(bins_r$wavg)

  # --- polynomial fit + pointwise CI bands -------------------------------------
  poly_fit <- function(xv, coefs, c0) {
    xc  <- xv - c0
    mat <- outer(xc, 0:length(coefs[-1]), "^")
    as.numeric(mat %*% coefs)
  }

  poly_se <- function(xv, V, c0) {
    xc   <- xv - c0
    mat  <- outer(xc, 0:(ncol(V) - 1), "^")
    vars <- rowSums((mat %*% V) * mat)
    sqrt(pmax(vars, 0))
  }

  xseq_l <- seq(c_val - h_l, c_val,       length.out = 200)
  xseq_r <- seq(c_val,       c_val + h_r, length.out = 200)
  yhat_l <- poly_fit(xseq_l, coefs_l, c_val)
  yhat_r <- poly_fit(xseq_r, coefs_r, c_val)
  fit_l  <- data.frame(x = xseq_l, y = yhat_l)
  fit_r  <- data.frame(x = xseq_r, y = yhat_r)

  if (ci) {
    se_l   <- poly_se(xseq_l, V_l, c_val)
    se_r   <- poly_se(xseq_r, V_r, c_val)
    band_l <- data.frame(x   = xseq_l,
                         ylo = yhat_l - z_crit * se_l,
                         yhi = yhat_l + z_crit * se_l)
    band_r <- data.frame(x   = xseq_r,
                         ylo = yhat_r - z_crit * se_r,
                         yhi = yhat_r + z_crit * se_r)
  }

  y0_l <- coefs_l[1]
  y0_r <- coefs_r[1]

  # --- RD plot annotation ------------------------------------------------------
  lbl <- sprintf("RD = %.3f\n%d%% RBC CI: [%.3f, %.3f]%s",
                 round(tau_cl,     3), lev,
                 round(ci_rb[[1]], 3),
                 round(ci_rb[[2]], 3),
                 if (is_fuzzy) "\n(Fuzzy)" else "")

  # --- build RD plot -----------------------------------------------------------
  p_rd <- ggplot2::ggplot()

  if (ci) {
    p_rd <- p_rd +
      ggplot2::geom_ribbon(data = band_l,
                           ggplot2::aes(x = x, ymin = ylo, ymax = yhi),
                           fill = col.l, alpha = 0.12) +
      ggplot2::geom_ribbon(data = band_r,
                           ggplot2::aes(x = x, ymin = ylo, ymax = yhi),
                           fill = col.r, alpha = 0.12)
  }

  p_rd <- p_rd +
    ggplot2::geom_point(data = bins_l,
                        ggplot2::aes(x = xmid, y = ymid, size = sz),
                        color = col.l, alpha = 0.85, show.legend = FALSE) +
    ggplot2::geom_point(data = bins_r,
                        ggplot2::aes(x = xmid, y = ymid, size = sz),
                        color = col.r, alpha = 0.85, show.legend = FALSE) +
    ggplot2::geom_line(data = fit_l,
                       ggplot2::aes(x = x, y = y),
                       color = col.l, linewidth = 1) +
    ggplot2::geom_line(data = fit_r,
                       ggplot2::aes(x = x, y = y),
                       color = col.r, linewidth = 1) +
    ggplot2::geom_vline(xintercept = c_val,
                        linetype = "dashed", color = "grey50", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(x = c_val, y = y0_l),
                        shape = 21, size = 3.5,
                        color = col.l, fill = "white", stroke = 1.5) +
    ggplot2::geom_point(ggplot2::aes(x = c_val, y = y0_r),
                        shape = 21, size = 3.5,
                        color = col.r, fill = "white", stroke = 1.5) +
    ggplot2::annotate("text",
                      x = c_val + 0.04 * (h_l + h_r), y = Inf,
                      label = lbl,
                      hjust = 0, vjust = 1.4,
                      size = txt_main, color = "grey20",
                      lineheight = 1.3) +
    ggplot2::scale_size_identity() +
    ggplot2::labs(title = title, x = x.label, y = y.label) +
    theme_rd()

  if (!is.null(x.lim)) p_rd <- p_rd + ggplot2::xlim(x.lim)
  if (!is.null(y.lim)) p_rd <- p_rd + ggplot2::ylim(y.lim)

  # --- effect plot: conventional estimate + RBC CI with significance stars -----
  if (show_effect) {
    pv_rb <- obj$pv["Robust", 1]
    stars <- ifelse(pv_rb < 0.01, "***",
             ifelse(pv_rb < 0.05, "**",
             ifelse(pv_rb < 0.10, "*", "")))

    eff_lbl <- sprintf("%.3f%s  [%.3f, %.3f]",
                       tau_cl, stars, ci_rb[[1]], ci_rb[[2]])

    eff_df <- data.frame(label = "RD Effect",
                         est   = tau_cl,
                         lo    = ci_rb[[1]],
                         hi    = ci_rb[[2]])

    eff_size <- base_size - 4
    eff_txt  <- eff_size * 0.353

    p_eff <- ggplot2::ggplot(eff_df,
                             ggplot2::aes(x = est, y = label)) +
      ggplot2::geom_vline(xintercept = 0,
                          linetype = "dashed", color = "grey60",
                          linewidth = 0.4) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = lo, xmax = hi),
                             width = 0, linewidth = 0.8, color = "grey30",
                             orientation = "y") +
      ggplot2::geom_point(size = 3, color = "grey20") +
      ggplot2::annotate("text",
                        x = tau_cl, y = 1.35,
                        label = eff_lbl,
                        size  = eff_txt,
                        color = "grey20") +
      ggplot2::labs(
        x       = sprintf("Estimate (%d%% RBC CI)", lev),
        y       = NULL,
        caption = "Point estimate: conventional. CI: robust bias-corrected. * p<0.10  ** p<0.05  *** p<0.01"
      ) +
      ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = c(0.4, 0.6))) +
      ggplot2::theme_bw(base_size = eff_size) +
      ggplot2::theme(
        panel.border        = ggplot2::element_blank(),
        axis.line.x         = ggplot2::element_line(color = "grey30",
                                                     linewidth = 0.4),
        panel.grid.major.x  = ggplot2::element_line(color = "grey92",
                                                     linewidth = 0.3),
        panel.grid.minor    = ggplot2::element_blank(),
        panel.grid.major.y  = ggplot2::element_blank(),
        axis.text.y         = ggplot2::element_text(size = eff_size,
                                                     color = "grey20"),
        axis.text.x         = ggplot2::element_text(size = eff_size - 1,
                                                     color = "grey30"),
        axis.title.x        = ggplot2::element_text(size = eff_size),
        axis.line.y         = ggplot2::element_blank(),
        axis.ticks.y        = ggplot2::element_blank(),
        plot.caption        = ggplot2::element_text(size = eff_size - 2,
                                                     color = "grey50",
                                                     hjust = 0),
        plot.margin         = ggplot2::margin(4, 16, 4, 8)
      )

    if (requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(p_rd, p_eff, ncol = 1, heights = c(3, 1))
    } else {
      print(p_rd)
      print(p_eff)
    }

    return(invisible(list(rd_plot = p_rd, effect_plot = p_eff)))
  }

  print(p_rd)
  invisible(p_rd)
}
