add_common_aes <- function(gplot, txtsize, scale_name = waiver(),
                           col = c("none", "full", "bw"),
                           col_aes = c("fill", "color"),
                           lval = 50,
                           greystart = 0.2,
                           greyend = 0.8,
                           continuous = c("none", "x", "y"),
                           n_x_ticks = 6,
                           n_y_ticks = 6,
                           xbreaks = NULL,
                           ybreaks = NULL,
                           xlim = NULL,
                           ylim = NULL,
                           xtrans = "identity",
                           ytrans = "identity",
                           xexpand = waiver(),
                           yexpand = waiver(),
                           facet_lab_txtsize = NULL,
                           ...) {
  p <- gplot +
    theme_bw() +
    theme(legend.title = element_text(size = txtsize),
          legend.text = element_text(size = txtsize - 3),
          title = element_text(face = "bold", size = (txtsize + 2)),
          axis.title.x = element_text(face = "bold", size = txtsize - 1),
          axis.title.y = element_text(face = "bold", size = txtsize - 1),
          axis.text.y = element_text(size = txtsize - 2),
          axis.text.x = element_text(size = txtsize - 2),
          strip.text.x = element_text(size = facet_lab_txtsize),
          strip.text.y = element_text(size = facet_lab_txtsize))

  col <- match.arg(col)
  col_aes <- match.arg(col_aes, several.ok = TRUE)
  if (col == "full") {
    if ("color" %in% col_aes) {
      p <- p +
        scale_color_discrete(name = scale_name, l = lval,
                             aesthetics = "color",
                             drop = FALSE)
    }
    if ("fill" %in% col_aes) {
      p <- p +
        scale_fill_discrete(name = scale_name, l = lval,
                            aesthetics = "fill",
                            drop = FALSE)
    }
  }
  if (col == "bw") {
    if ("color" %in% col_aes) {
      p <- p +
        scale_color_grey(name = scale_name, start = greystart, end = greyend,
                         aesthetics = "color",
                         drop = FALSE)
    }
    if ("fill" %in% col_aes) {
      p <- p +
        scale_fill_grey(name = scale_name, start = greystart, end = greyend,
                        aesthetics = "fill",
                        drop = FALSE)
    }
  }

  # axes and axis ticks
  continuous <- match.arg(continuous, several.ok = TRUE)

  if ("x" %in% continuous) {
    if (!is.null(xbreaks)) {
      xb <- xbreaks
    } else {
      xb <- number_ticks(n_x_ticks)
    }
    p <- p +
      scale_x_continuous(breaks = xb,
                         labels = scales::label_comma(),
                         limits = xlim,
                         trans = xtrans,
                         expand = xexpand)
  }
  if ("y" %in% continuous) {
    if (!is.null(ybreaks)) {
      yb <- ybreaks
    } else {
      yb <- number_ticks(n_y_ticks)
    }
    p <- p +
      scale_y_continuous(breaks = yb,
                         labels = labfun,
                         limits = ylim,
                         trans = ytrans,
                         expand = yexpand)
  }
  return(p)
}

library(formattable)
#' used to automatically label continuous scales
#' @keywords internal
#' @param x axis breaks
#' @return  a character vector giving a label for each input value
labfun <- function(x) {
  if (any(x > 999, na.rm = TRUE)) {
    comma(x, digits = 0)
  } else {
    x
  }
}


#' Number of ticks for \code{ggplot2} plots
#'
#' Function for determining number of ticks on axis of \code{ggplot2} plots.
#' @param n integer giving the desired number of ticks on axis of
#' \code{ggplot2} plots. Non-integer values are rounded down.
#' @section Details:
#' Based on function \code{pretty}.
#' @return a vector of axis-label breaks
#' @export
number_ticks <- function(n) {
  function(limits) {
    pretty(limits, n + 1)
  }
}

plot_icers_mod <- function(x,
                           txtsize = 12,
                           currency = "$",
                           effect_units = "QALYs",
                           label = c("frontier", "all", "none"),
                           label_max_char = NULL,
                           plot_frontier_only = FALSE,
                           alpha = 1,
                           n_x_ticks = 6,
                           n_y_ticks = 6,
                           xbreaks = NULL,
                           ybreaks = NULL,
                           xlim = NULL,
                           ylim = NULL,
                           xexpand = expansion(0.1),
                           yexpand = expansion(0.1),
                           max.iter = 20000,
                           ...) {
  if (ncol(x) > 7) {
    # reformat icers class object if uncertainty bounds are present
    x <- x %>%
      select(.data$Strategy, .data$Cost, .data$Effect,
             .data$Inc_Cost, .data$Inc_Effect,
             .data$ICER, .data$Status)
  }

  # type checking
  label <- match.arg(label)

  # this is so non-dominated strategies are plotted last (on top)
  x <- arrange(x, .data$Status)

  # change status text in data frame for plotting
  d_name <- "Dominated"
  ed_name <- "Weakly Dominated"
  nd_name <- "Efficient Frontier"

  status_expand <- c("D" = d_name, "ED" = ed_name,
                     "ND" = nd_name, "ref" = nd_name)
  x$Status <- factor(status_expand[x$Status], ordered = FALSE,
                     levels = c(d_name, ed_name, nd_name))

  # linetype
  plot_lines <- c("Dominated" = "blank",
                  "Weakly Dominated" = "blank",
                  "Efficient Frontier" = "solid")

  # names to refer to in aes_
  stat_name <- "Status"
  strat_name <- "Strategy"
  eff_name <- "Effect"
  cost_name <- "Cost"

  # frontier only
  if (plot_frontier_only) {
    plt_data <- x[x$Status == nd_name, ]
  } else {
    plt_data <- x
  }

  # make plot
  icer_plot <- ggplot(plt_data, aes_(x = as.name(eff_name),
                                     y = as.name(cost_name),
                                     shape = as.name(stat_name))) +
    geom_point(aes_(color = as.name(strat_name)), alpha = alpha, size = 2) +
    geom_line(aes_(linetype = as.name(stat_name), group = as.name(stat_name))) +
    scale_linetype_manual(name = NULL, values = plot_lines) +
    scale_shape_discrete(name = NULL) +
    scale_color_nejm() +
    labs(x = paste0("Effect (", effect_units, ")"),
         y = paste0("Cost (", currency, ")"))

  icer_plot <- add_common_aes(icer_plot, txtsize, col = "none",
                              continuous = c("x", "y"),
                              n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                              xbreaks = xbreaks, ybreaks = ybreaks,
                              xlim = xlim, ylim = ylim,
                              xexpand = xexpand, yexpand = yexpand)

  # labeling
  if (label != "none") {
    if (!is.null(label_max_char)) {
      plt_data[, strat_name] <- str_sub(plt_data[, strat_name],
                                        start = 1L, end = label_max_char)
    }
    if (label == "all") {
      lab_data <- plt_data
    }
    if (label == "frontier") {
      lab_data <- plt_data[plt_data$Status == nd_name, ]
    }

    icer_plot <- icer_plot +
      geom_label_repel(data = lab_data,
                       aes_(label = as.name(strat_name)),
                       size = 3,
                       show.legend = FALSE,
                       max.iter = max.iter,
                       direction = "both")
  }
  return(icer_plot)
}

plot_icers_mod_subgrp <- function(x,
                                  txtsize = 12,
                                  currency = "$",
                                  effect_units = "QALYs",
                                  label = c("frontier", "all", "none"),
                                  label_max_char = NULL,
                                  plot_frontier_only = FALSE,
                                  alpha = 1,
                                  n_x_ticks = 6,
                                  n_y_ticks = 6,
                                  xbreaks = NULL,
                                  ybreaks = NULL,
                                  xlim = NULL,
                                  ylim = NULL,
                                  xexpand = expansion(0.1),
                                  yexpand = expansion(0.1),
                                  max.iter = 20000,
                                  ...) {
  # if (ncol(x) > 7) {
  #   # reformat icers class object if uncertainty bounds are present
  #   x <- x %>%
  #     select(.data$Strategy, .data$Cost, .data$Effect,
  #            .data$Inc_Cost, .data$Inc_Effect,
  #            .data$ICER, .data$Status)
  # }

  # type checking
  label <- match.arg(label)

  # this is so non-dominated strategies are plotted last (on top)
  x <- arrange(x, .data$Status)

  # change status text in data frame for plotting
  d_name <- "Dominated"
  ed_name <- ""
  nd_name <- "Efficient Frontier"

  status_expand <- c("D" = d_name, "ED" = ed_name,
                     "ND" = nd_name, "ref" = nd_name)
  x$Status <- factor(status_expand[x$Status], ordered = FALSE,
                     levels = c(d_name, ed_name, nd_name))

  # linetype
  plot_lines <- c("Dominated" = "blank",
                  "Weakly Dominated" = "blank",
                  "Efficient Frontier" = "solid")

  # names to refer to in aes_
  stat_name <- "Status"
  strat_name <- "Strategy"
  eff_name <- "Effect"
  cost_name <- "Cost"
  sub_name <- "subgrp"

  # frontier only
  if (plot_frontier_only) {
    plt_data <- x[x$Status == nd_name, ]
  } else {
    plt_data <- x
  }
  l_plt_data <- split(plt_data, plt_data$subgrp)

  print(l_plt_data[[1]])

  # make plot
  icer_plot <- ggplot(l_plt_data[[1]], aes_(x = as.name(eff_name),
                                           y = as.name(cost_name))) +
    geom_point(aes_(shape = as.name(strat_name),
                   color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(aes_(linetype = as.name(stat_name), group = as.name(stat_name))) +
    geom_point(data = l_plt_data[[2]], aes_(shape = as.name(strat_name),
                                            color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(data = l_plt_data[[2]], aes_(x = as.name(eff_name),
                                          y = as.name(cost_name),
                                          linetype = as.name(stat_name), group = as.name(stat_name))) +
    geom_point(data = l_plt_data[[3]], aes_(x = as.name(eff_name),
                                           y = as.name(cost_name),
                                           shape = as.name(strat_name),
                                           color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(data = l_plt_data[[3]], aes_(x = as.name(eff_name),
                                          y = as.name(cost_name),
                                          linetype = as.name(stat_name), group = as.name(stat_name))) +
    geom_point(data = l_plt_data[[4]], aes_(x = as.name(eff_name),
                                            y = as.name(cost_name),
                                            shape = as.name(strat_name),
                                            color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(data = l_plt_data[[4]], aes_(x = as.name(eff_name),
                                           y = as.name(cost_name),
                                           linetype = as.name(stat_name), group = as.name(stat_name))) +
    scale_linetype_manual(name = NULL, values = plot_lines) +
    scale_shape_manual(name = NULL, values = c(16, 17, 15, 18, 3, 8)) +
    scale_color_nejm() +
    labs(x = paste0("Effect (", effect_units, ")"),
         y = paste0("Cost (", currency, ")"))

  icer_plot <- add_common_aes(icer_plot, txtsize, col = "none",
                              continuous = c("x", "y"),
                              n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                              xbreaks = xbreaks, ybreaks = ybreaks,
                              xlim = xlim, ylim = ylim,
                              xexpand = xexpand, yexpand = yexpand)

  # labeling
  if (label != "none") {
    if (!is.null(label_max_char)) {
      plt_data[, strat_name] <- str_sub(plt_data[, strat_name],
                                        start = 1L, end = label_max_char)
    }
    if (label == "all") {
      lab_data <- plt_data
    }
    if (label == "frontier") {
      lab_data <- plt_data[plt_data$Status == nd_name, ]
    }

    icer_plot <- icer_plot +
      geom_label_repel(data = lab_data,
                       aes_(label = as.name(strat_name)),
                       size = 3,
                       show.legend = FALSE,
                       max.iter = max.iter,
                       direction = "both")
  }
  return(icer_plot)
}

plot.ceac <- function(x,
                      frontier = TRUE,
                      points = TRUE,
                      currency = "$",
                      min_prob = 0,
                      txtsize = 12,
                      n_x_ticks = 10,
                      n_y_ticks = 8,
                      xbreaks = NULL,
                      ybreaks = NULL,
                      ylim = NULL,
                      xlim = c(0, NA),
                      col = c("full", "bw"),
                      ...) {
  wtp_name <- "WTP"
  prop_name <- "Proportion"
  strat_name <- "Strategy"
  x$WTP_thou <- x[, wtp_name] / 1000

  # removing strategies with probabilities always below `min_prob`
  # get group-wise max probability
  if (min_prob > 0) {
    max_prob <- x %>%
      group_by(.data$Strategy) %>%
      summarize(maxpr = max(.data$Proportion)) %>%
      filter(.data$maxpr >= min_prob)
    strat_to_keep <- max_prob$Strategy
    if (length(strat_to_keep) == 0) {
      stop(
        paste("no strategies remaining. you may want to lower your min_prob value (currently ",
              min_prob, ")", sep = "")
      )
    }
    # report filtered out strategies
    old_strat <- unique(x$Strategy)
    diff_strat <- setdiff(old_strat, strat_to_keep)
    n_diff_strat <- length(diff_strat)
    if (n_diff_strat > 0) {
      # report strategies filtered out
      cat("filtered out ", n_diff_strat, " strategies with max prob below ", min_prob, ":\n",
          paste(diff_strat, collapse = ","), "\n", sep = "")

      # report if any filtered strategies are on the frontier
      df_filt <- filter(x, .data$Strategy %in% diff_strat & .data$On_Frontier)
      if (nrow(df_filt) > 0) {
        cat(paste0("WARNING - some strategies that were filtered out are on the frontier:\n",
                   paste(unique(df_filt$Strategy), collapse = ","), "\n"))
      }
    }

    # filter dataframe
    x <- filter(x, .data$Strategy %in% strat_to_keep)
  }

  # Drop unused strategy names
  x$Strategy <- droplevels(x$Strategy)

  p <- ggplot(data = x, aes_(x = as.name("WTP_thou"),
                             y = as.name(prop_name),
                             color = as.name(strat_name))) +
    geom_line() +
    theme_bw() +
    xlab(paste("Willingness to Pay (Thousand ", currency, " / QALY)", sep = "")) +
    ylab("Pr Cost-Effective")

  p <- p + geom_point(aes_(color = as.name(strat_name)))

  front <- x[x$On_Frontier, ]
  p <- p + geom_point(data = front, aes_(x = as.name("WTP_thou"),
                                         y = as.name(prop_name),
                                         shape = as.name("On_Frontier")),
                      size = 3, stroke = 1, color = "black") +
    scale_shape_manual(name = NULL, values = 0, labels = "Frontier") +
    guides(color = guide_legend(order = 1),
           shape = guide_legend(order = 2))

  col <- match.arg(col)
  add_common_aes(p, txtsize, col = col, col_aes = "color",
                 continuous = c("x", "y"), n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                 xbreaks = xbreaks, ybreaks = ybreaks,
                 ylim = ylim, xlim = xlim)
  p <- p + scale_color_nejm(labels = c("Status Quo", "5% Increase",
                                       "10% Increase", "15% Increase",
                                       "20% Increase", "25% Increase"))
  return(p)
}

plot.exp_loss <- function(x,
                          log_y = TRUE,
                          frontier = TRUE,
                          points = TRUE,
                          lsize = 1,
                          txtsize = 12,
                          currency = "$",
                          effect_units = "QALY",
                          n_y_ticks = 8,
                          n_x_ticks = 20,
                          xbreaks = NULL,
                          ybreaks = NULL,
                          xlim = c(0, NA),
                          ylim = NULL,
                          col = c("full", "bw"),
                          ...) {
  wtp_name <- "WTP_thou"
  loss_name <- "Expected_Loss"
  strat_name <- "Strategy"
  x[, wtp_name] <- x$WTP / 1000

  # split into on frontier and not on frontier
  nofront <- x
  front <- x[x$On_Frontier, ]

  # Drop unused levels from strategy names
  nofront$Strategy <- droplevels(nofront$Strategy)
  front$Strategy <- droplevels(front$Strategy)
  # formatting if logging the y axis
  if (log_y) {
    tr <- "log10"
  } else {
    tr <- "identity"
  }

  p <- ggplot(data = nofront, aes_(x = as.name(wtp_name),
                                   y = as.name(loss_name))) +
    xlab(paste0("Willingness to Pay (Thousand ", currency, "/", effect_units, ")")) +
    ylab(paste0("Expected Loss (", currency, ")"))

  # color
  col <- match.arg(col)
  ## change linetype too if color is black and white
  if (col == "full") {
    if (points) {
      p <- p + geom_point(aes_(color = as.name(strat_name)))
    }
    p <- p +
      geom_line(size = lsize, aes_(color = as.name(strat_name)))

  }
  if (col == "bw") {
    if (points) {
      p <- p + geom_point()
    }
    p <- p +
      geom_line(aes_(linetype = as.name(strat_name)))
  }

  p <- add_common_aes(p, txtsize, col = col, col_aes = c("color", "line"),
                      continuous = c("x", "y"),
                      n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                      xbreaks = xbreaks, ybreaks = ybreaks,
                      xlim = xlim, ylim = ylim,
                      ytrans = tr)
  if (frontier) {
    p <- p + geom_point(data = front, aes_(x = as.name(wtp_name),
                                           y = as.name(loss_name),
                                           shape = as.name("On_Frontier")),
                        size = 3, stroke = 1, color = "black") +
      scale_shape_manual(name = NULL, values = 0, labels = "Frontier") +
      guides(color = guide_legend(order = 1),
             linetype = guide_legend(order = 1),
             shape = guide_legend(order = 2))
  }
  p <- p + scale_color_nejm(labels = c("Status Quo", "5% Increase",
                                       "10% Increase", "15% Increase",
                                       "20% Increase", "25% Increase"))
  return(p)
}

plot_icers_mod_subgrp_2 <- function(x,
                                  txtsize = 12,
                                  currency = "$",
                                  effect_units = "QALYs",
                                  label = c("frontier", "all", "none"),
                                  label_max_char = NULL,
                                  plot_frontier_only = FALSE,
                                  alpha = 1,
                                  n_x_ticks = 6,
                                  n_y_ticks = 6,
                                  xbreaks = NULL,
                                  ybreaks = NULL,
                                  xlim = NULL,
                                  ylim = NULL,
                                  xexpand = expansion(0.1),
                                  yexpand = expansion(0.1),
                                  max.iter = 20000,
                                  ...) {
  # if (ncol(x) > 7) {
  #   # reformat icers class object if uncertainty bounds are present
  #   x <- x %>%
  #     select(.data$Strategy, .data$Cost, .data$Effect,
  #            .data$Inc_Cost, .data$Inc_Effect,
  #            .data$ICER, .data$Status)
  # }

  # type checking
  label <- match.arg(label)

  # this is so non-dominated strategies are plotted last (on top)
  x <- arrange(x, .data$Status)

  # change status text in data frame for plotting
  d_name <- "Dominated"
  ed_name <- ""
  nd_name <- "Efficient Frontier"

  status_expand <- c("D" = d_name, "ED" = ed_name,
                     "ND" = nd_name, "ref" = nd_name)
  x$Status <- factor(status_expand[x$Status], ordered = FALSE,
                     levels = c(d_name, ed_name, nd_name))

  # linetype
  plot_lines <- c("Dominated" = "blank",
                  "Weakly Dominated" = "blank",
                  "Efficient Frontier" = "solid")

  # names to refer to in aes_
  stat_name <- "Status"
  strat_name <- "Strategy"
  eff_name <- "Effect"
  cost_name <- "Cost"
  sub_name <- "subgrp"

  # frontier only
  if (plot_frontier_only) {
    plt_data <- x[x$Status == nd_name, ]
  } else {
    plt_data <- x
  }
  l_plt_data <- split(plt_data, plt_data$subgrp)

  print(l_plt_data[[1]])

  # make plot
  icer_plot <- ggplot(l_plt_data[[1]], aes_(x = as.name(eff_name),
                                            y = as.name(cost_name))) +
    geom_point(aes_(shape = as.name(strat_name),
                    color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(aes_(linetype = as.name(stat_name), group = as.name(stat_name))) +
    geom_point(data = l_plt_data[[2]], aes_(shape = as.name(strat_name),
                                            color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(data = l_plt_data[[2]], aes_(x = as.name(eff_name),
                                           y = as.name(cost_name),
                                           linetype = as.name(stat_name), group = as.name(stat_name))) +
    scale_linetype_manual(name = NULL, values = plot_lines) +
    scale_shape_manual(name = NULL, values = c(16, 17, 15, 18, 3, 8)) +
    scale_color_nejm() +
    labs(x = paste0("Effect (", effect_units, ")"),
         y = paste0("Cost (", currency, ")"))

  icer_plot <- add_common_aes(icer_plot, txtsize, col = "none",
                              continuous = c("x", "y"),
                              n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                              xbreaks = xbreaks, ybreaks = ybreaks,
                              xlim = xlim, ylim = ylim,
                              xexpand = xexpand, yexpand = yexpand)

  # labeling
  if (label != "none") {
    if (!is.null(label_max_char)) {
      plt_data[, strat_name] <- str_sub(plt_data[, strat_name],
                                        start = 1L, end = label_max_char)
    }
    if (label == "all") {
      lab_data <- plt_data
    }
    if (label == "frontier") {
      lab_data <- plt_data[plt_data$Status == nd_name, ]
    }

    icer_plot <- icer_plot +
      geom_label_repel(data = lab_data,
                       aes_(label = as.name(strat_name)),
                       size = 3,
                       show.legend = FALSE,
                       max.iter = max.iter,
                       direction = "both")
  }
  return(icer_plot)
}

plot_icers_mod_subgrp_3 <- function(x,
                                  txtsize = 12,
                                  currency = "$",
                                  effect_units = "QALYs",
                                  label = c("frontier", "all", "none"),
                                  label_max_char = NULL,
                                  plot_frontier_only = FALSE,
                                  alpha = 1,
                                  n_x_ticks = 6,
                                  n_y_ticks = 6,
                                  xbreaks = NULL,
                                  ybreaks = NULL,
                                  xlim = NULL,
                                  ylim = NULL,
                                  xexpand = expansion(0.1),
                                  yexpand = expansion(0.1),
                                  max.iter = 20000,
                                  ...) {
  # if (ncol(x) > 7) {
  #   # reformat icers class object if uncertainty bounds are present
  #   x <- x %>%
  #     select(.data$Strategy, .data$Cost, .data$Effect,
  #            .data$Inc_Cost, .data$Inc_Effect,
  #            .data$ICER, .data$Status)
  # }

  # type checking
  label <- match.arg(label)

  # this is so non-dominated strategies are plotted last (on top)
  x <- arrange(x, .data$Status)

  # change status text in data frame for plotting
  d_name <- "Dominated"
  ed_name <- ""
  nd_name <- "Efficient Frontier"

  status_expand <- c("D" = d_name, "ED" = ed_name,
                     "ND" = nd_name, "ref" = nd_name)
  x$Status <- factor(status_expand[x$Status], ordered = FALSE,
                     levels = c(d_name, ed_name, nd_name))

  # linetype
  plot_lines <- c("Dominated" = "blank",
                  "Weakly Dominated" = "blank",
                  "Efficient Frontier" = "solid")

  # names to refer to in aes_
  stat_name <- "Status"
  strat_name <- "Strategy"
  eff_name <- "Effect"
  cost_name <- "Cost"
  sub_name <- "subgrp"

  # frontier only
  if (plot_frontier_only) {
    plt_data <- x[x$Status == nd_name, ]
  } else {
    plt_data <- x
  }
  l_plt_data <- split(plt_data, plt_data$subgrp)

  print(l_plt_data[[1]])

  # make plot
  icer_plot <- ggplot(l_plt_data[[1]], aes_(x = as.name(eff_name),
                                            y = as.name(cost_name))) +
    geom_point(aes_(shape = as.name(strat_name),
                    color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(aes_(linetype = as.name(stat_name), group = as.name(stat_name))) +
    geom_point(data = l_plt_data[[2]], aes_(shape = as.name(strat_name),
                                            color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(data = l_plt_data[[2]], aes_(x = as.name(eff_name),
                                           y = as.name(cost_name),
                                           linetype = as.name(stat_name), group = as.name(stat_name))) +
    geom_point(data = l_plt_data[[3]], aes_(x = as.name(eff_name),
                                            y = as.name(cost_name),
                                            shape = as.name(strat_name),
                                            color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(data = l_plt_data[[3]], aes_(x = as.name(eff_name),
                                           y = as.name(cost_name),
                                           linetype = as.name(stat_name), group = as.name(stat_name))) +
    scale_linetype_manual(name = NULL, values = plot_lines) +
    scale_shape_manual(name = NULL, values = c(16, 17, 15, 18, 3, 8)) +
    scale_color_nejm() +
    labs(x = paste0("Effect (", effect_units, ")"),
         y = paste0("Cost (", currency, ")"))

  icer_plot <- add_common_aes(icer_plot, txtsize, col = "none",
                              continuous = c("x", "y"),
                              n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                              xbreaks = xbreaks, ybreaks = ybreaks,
                              xlim = xlim, ylim = ylim,
                              xexpand = xexpand, yexpand = yexpand)

  # labeling
  if (label != "none") {
    if (!is.null(label_max_char)) {
      plt_data[, strat_name] <- str_sub(plt_data[, strat_name],
                                        start = 1L, end = label_max_char)
    }
    if (label == "all") {
      lab_data <- plt_data
    }
    if (label == "frontier") {
      lab_data <- plt_data[plt_data$Status == nd_name, ]
    }

    icer_plot <- icer_plot +
      geom_label_repel(data = lab_data,
                       aes_(label = as.name(strat_name)),
                       size = 3,
                       show.legend = FALSE,
                       max.iter = max.iter,
                       direction = "both")
  }
  return(icer_plot)
}

plot_icers_mod_subgrp_4 <- function(x,
                                  txtsize = 12,
                                  currency = "$",
                                  effect_units = "QALYs",
                                  label = c("frontier", "all", "none"),
                                  label_max_char = NULL,
                                  plot_frontier_only = FALSE,
                                  alpha = 1,
                                  n_x_ticks = 6,
                                  n_y_ticks = 6,
                                  xbreaks = NULL,
                                  ybreaks = NULL,
                                  xlim = NULL,
                                  ylim = NULL,
                                  xexpand = expansion(0.1),
                                  yexpand = expansion(0.1),
                                  max.iter = 20000,
                                  ...) {
  # if (ncol(x) > 7) {
  #   # reformat icers class object if uncertainty bounds are present
  #   x <- x %>%
  #     select(.data$Strategy, .data$Cost, .data$Effect,
  #            .data$Inc_Cost, .data$Inc_Effect,
  #            .data$ICER, .data$Status)
  # }

  # type checking
  label <- match.arg(label)

  # this is so non-dominated strategies are plotted last (on top)
  x <- arrange(x, .data$Status)

  # change status text in data frame for plotting
  d_name <- "Dominated"
  ed_name <- ""
  nd_name <- "Efficient Frontier"

  status_expand <- c("D" = d_name, "ED" = ed_name,
                     "ND" = nd_name, "ref" = nd_name)
  x$Status <- factor(status_expand[x$Status], ordered = FALSE,
                     levels = c(d_name, ed_name, nd_name))

  # linetype
  plot_lines <- c("Dominated" = "blank",
                  "Weakly Dominated" = "blank",
                  "Efficient Frontier" = "solid")

  # names to refer to in aes_
  stat_name <- "Status"
  strat_name <- "Strategy"
  eff_name <- "Effect"
  cost_name <- "Cost"
  sub_name <- "subgrp"

  # frontier only
  if (plot_frontier_only) {
    plt_data <- x[x$Status == nd_name, ]
  } else {
    plt_data <- x
  }
  l_plt_data <- split(plt_data, plt_data$subgrp)

  print(l_plt_data[[1]])

  # make plot
  icer_plot <- ggplot(l_plt_data[[1]], aes_(x = as.name(eff_name),
                                            y = as.name(cost_name))) +
    geom_point(aes_(shape = as.name(strat_name),
                    color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(aes_(linetype = as.name(stat_name), group = as.name(stat_name))) +
    geom_point(data = l_plt_data[[2]], aes_(shape = as.name(strat_name),
                                            color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(data = l_plt_data[[2]], aes_(x = as.name(eff_name),
                                           y = as.name(cost_name),
                                           linetype = as.name(stat_name), group = as.name(stat_name))) +
    geom_point(data = l_plt_data[[3]], aes_(x = as.name(eff_name),
                                            y = as.name(cost_name),
                                            shape = as.name(strat_name),
                                            color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(data = l_plt_data[[3]], aes_(x = as.name(eff_name),
                                           y = as.name(cost_name),
                                           linetype = as.name(stat_name), group = as.name(stat_name))) +
    geom_point(data = l_plt_data[[4]], aes_(x = as.name(eff_name),
                                            y = as.name(cost_name),
                                            shape = as.name(strat_name),
                                            color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(data = l_plt_data[[4]], aes_(x = as.name(eff_name),
                                           y = as.name(cost_name),
                                           linetype = as.name(stat_name), group = as.name(stat_name))) +
    scale_linetype_manual(name = NULL, values = plot_lines) +
    scale_shape_manual(name = NULL, values = c(16, 17, 15, 18, 3, 8)) +
    scale_color_nejm() +
    labs(x = paste0("Effect (", effect_units, ")"),
         y = paste0("Cost (", currency, ")"))

  icer_plot <- add_common_aes(icer_plot, txtsize, col = "none",
                              continuous = c("x", "y"),
                              n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                              xbreaks = xbreaks, ybreaks = ybreaks,
                              xlim = xlim, ylim = ylim,
                              xexpand = xexpand, yexpand = yexpand)

  # labeling
  if (label != "none") {
    if (!is.null(label_max_char)) {
      plt_data[, strat_name] <- str_sub(plt_data[, strat_name],
                                        start = 1L, end = label_max_char)
    }
    if (label == "all") {
      lab_data <- plt_data
    }
    if (label == "frontier") {
      lab_data <- plt_data[plt_data$Status == nd_name, ]
    }

    icer_plot <- icer_plot +
      geom_label_repel(data = lab_data,
                       aes_(label = as.name(strat_name)),
                       size = 3,
                       show.legend = FALSE,
                       max.iter = max.iter,
                       direction = "both")
  }
  return(icer_plot)
}

plot_icers_mod_subgrp <- function(x,
                                  txtsize = 12,
                                  currency = "$",
                                  effect_units = "QALYs",
                                  label = c("frontier", "all", "none"),
                                  label_max_char = NULL,
                                  plot_frontier_only = FALSE,
                                  alpha = 1,
                                  n_x_ticks = 6,
                                  n_y_ticks = 6,
                                  xbreaks = NULL,
                                  ybreaks = NULL,
                                  xlim = NULL,
                                  ylim = NULL,
                                  xexpand = expansion(0.1),
                                  yexpand = expansion(0.1),
                                  max.iter = 20000,
                                  ...) {
  # if (ncol(x) > 7) {
  #   # reformat icers class object if uncertainty bounds are present
  #   x <- x %>%
  #     select(.data$Strategy, .data$Cost, .data$Effect,
  #            .data$Inc_Cost, .data$Inc_Effect,
  #            .data$ICER, .data$Status)
  # }

  # type checking
  label <- match.arg(label)

  # this is so non-dominated strategies are plotted last (on top)
  x <- arrange(x, .data$Status)

  # change status text in data frame for plotting
  d_name <- "Dominated"
  ed_name <- ""
  nd_name <- "Efficient Frontier"

  status_expand <- c("D" = d_name, "ED" = ed_name,
                     "ND" = nd_name, "ref" = nd_name)
  x$Status <- factor(status_expand[x$Status], ordered = FALSE,
                     levels = c(d_name, ed_name, nd_name))

  # linetype
  plot_lines <- c("Dominated" = "blank",
                  "Weakly Dominated" = "blank",
                  "Efficient Frontier" = "solid")
  plot_lines <- c("Efficient Frontier" = "solid")
  # names to refer to in aes_
  stat_name <- "Status"
  strat_name <- "Strategy"
  eff_name <- "Effect"
  cost_name <- "Cost"
  sub_name <- "subgrp"

  # frontier only
  if (plot_frontier_only) {
    plt_data <- x[x$Status == nd_name, ]
  } else {
    plt_data <- x
  }
  l_plt_data <- split(plt_data, plt_data$subgrp)

  # plt_data$subgrp <- factor(plt_data$subgrp,
  #                           levels = c("5%", "10%", "15%", "20%"))

  # make plot
  # icer_plot <- ggplot(l_plt_data[[1]], aes_(x = as.name(eff_name),
  #                                           y = as.name(cost_name))) +
  #   geom_point(aes_(shape = as.name(strat_name),
  #                   color = as.name(sub_name)), alpha = alpha, size = 2) +
  #   geom_line(aes_(linetype = as.name(stat_name), group = as.name(stat_name))) +
  #   geom_point(data = l_plt_data[[2]], aes_(shape = as.name(strat_name),
  #                                           color = as.name(sub_name)), alpha = alpha, size = 2) +
  #   geom_line(data = l_plt_data[[2]], aes_(x = as.name(eff_name),
  #                                          y = as.name(cost_name),
  #                                          linetype = as.name(stat_name), group = as.name(stat_name))) +
  #   geom_point(data = l_plt_data[[3]], aes_(x = as.name(eff_name),
  #                                           y = as.name(cost_name),
  #                                           shape = as.name(strat_name),
  #                                           color = as.name(sub_name)), alpha = alpha, size = 2) +
  #   geom_line(data = l_plt_data[[3]], aes_(x = as.name(eff_name),
  #                                          y = as.name(cost_name),
  #                                          linetype = as.name(stat_name), group = as.name(stat_name))) +
  #   geom_point(data = l_plt_data[[4]], aes_(x = as.name(eff_name),
  #                                           y = as.name(cost_name),
  #                                           shape = as.name(strat_name),
  #                                           color = as.name(sub_name)), alpha = alpha, size = 2) +
  #   geom_line(data = l_plt_data[[4]], aes_(x = as.name(eff_name),
  #                                          y = as.name(cost_name),
  #                                          linetype = as.name(stat_name), group = as.name(stat_name))) +
  #   scale_linetype_manual(name = NULL, values = plot_lines) +
  #   scale_shape_manual(name = NULL, values = c(16, 17, 15, 18, 3, 8)) +
  #   scale_color_nejm() +
  #   labs(x = paste0("Effect (", effect_units, ")"),
  #        y = paste0("Cost (", currency, ")"))

  icer_plot <- ggplot(plt_data, aes_(y = as.name(eff_name),
                                     x = as.name(cost_name))) +
    geom_point(aes_(shape = as.name(strat_name),
                    color = as.name(sub_name)), alpha = alpha, size = 2) +
    geom_line(aes_(linetype = as.name(stat_name), group = as.name(stat_name))) +
    facet_grid(as.factor(subgrp) ~ Perspective, scales = "free_x") +
    scale_linetype_manual(name = NULL, values = plot_lines) +
    scale_color_nejm() +
    labs(y = paste0("Effect (", effect_units, ")"),
         x = paste0("Cost (", currency, ")"))


  icer_plot <- add_common_aes(icer_plot, txtsize, col = "none",
                              continuous = c("x", "y"),
                              n_x_ticks = n_x_ticks, n_y_ticks = n_y_ticks,
                              xbreaks = xbreaks, ybreaks = ybreaks,
                              xlim = xlim, ylim = ylim,
                              xexpand = xexpand, yexpand = yexpand)
  # labeling
  if (label != "none") {
    if (!is.null(label_max_char)) {
      plt_data[, strat_name] <- str_sub(plt_data[, strat_name],
                                        start = 1L, end = label_max_char)
    }
    if (label == "all") {
      lab_data <- plt_data
    }
    if (label == "frontier") {
      lab_data <- plt_data[plt_data$Status == nd_name, ]
    }

    icer_plot <- icer_plot +
      geom_label_repel(data = lab_data,
                       aes_(label = as.name(strat_name)),
                       size = 3,
                       show.legend = FALSE,
                       max.iter = max.iter,
                       direction = "both")
  }
  return(icer_plot)
}
