#' Plot surveysd-Objects
#'
#' Plot results of `calc.stError()`
#'
#' @param x object of class 'surveysd' output of function [calc.stError]
#' @param variable Name of the variable for which standard errors have been
#'   calcualated in `dat`
#' @param type can bei either `"summary"` or `"grouping"`, default value is
#'   `"summary"`. For `"summary"` a barplot is created giving an overview of the
#'   number of estimates having the flag `smallGroup`, `cvHigh`, both or none
#'   of them. For 'grouping' results for point estimate and standard error are
#'   plotted for pre defined groups.
#' @param groups If `type='grouping'` variables must be defined by which the
#'   data is grouped. Only 2 levels are supported as of right now. If only one
#'   group is defined the higher group will be the estimate over the whole
#'   period. Results are plotted for the first argument in `groups` as well as
#'   for the combination of `groups[1]` and `groups[2]`.
#' @param sd.type can bei either `'ribbon'` or `'dot'` and is only used if
#'   `type='grouping'`. Default is `"dot"`
#' For `sd.type='dot'` point estimates are plotted and flagged if the
#' corresponding standard error and/or the standard error using the mean over
#' k-periods exceeded the value `cv.limit` (see [calc.stError]).
#' For `sd.type='ribbon'` the point estimates including ribbons, defined by
#' point estimate +- estimated standard error are plotted.
#' The calculated standard errors using the mean over k periods are plotted
#' using less transparency. Results for the higher level (~`groups[1]`) are
#' coloured grey.
#' @param ... additional arguments supplied to plot.
#'
#' @examples
#'
#' library(surveysd)
#' 
#' set.seed(1234)
#' eusilc <- demo.eusilc(n = 3, prettyNames = TRUE)
#'
#' dat_boot <- draw.bootstrap(eusilc, REP = 3, hid = "hid", weights = "pWeight",
#'                            strata = "region", period = "year")
#'
#' # calibrate weight for bootstrap replicates
#' dat_boot_calib <- recalib(dat_boot, conP.var = "gender", conH.var = "region")
#'
#' # estimate weightedRatio for povmd60 per period
#' group <- list("gender", "region", c("gender", "region"))
#' err.est <- calc.stError(dat_boot_calib, var = "povertyRisk",
#'                         fun = weightedRatio,
#'                         group = group , period.mean = NULL)
#'
#'
#' plot(err.est)
#'
#' # plot results for gender
#' # dotted line is the result on the national level
#' plot(err.est, type = "grouping", groups = "gender")
#'
#'
#' # plot results for rb090 in each db040
#' # with standard errors as ribbons
#' plot(err.est, type = "grouping", groups = c("gender", "region"), sd.type = "ribbon")
#'
#'
#' @export

plot.surveysd <- function(
  x, variable = x$param$var[1], type = c("summary", "grouping"), groups =
    NULL, sd.type = c("dot", "ribbon"), ...) {

  res_type <- GROUPING <- shape_bool <- shape_bool2 <- NULL


  #################
  # Input checking
  type <- type[1]
  sd.type <- sd.type[1]
  if (class(x)[1] != "surveysd") {
    stop("The data needs to be an object of class 'surveysd'!")
  }
  if (!type %in% c("summary", "grouping")) {
    stop("Parameter type can only take values 'summary' or 'grouping'!")
  }
  if (!variable %in% x$param$var) {
    stop("No results for ", variable, " present in the data!")
  }
  if (type == "grouping") {

    if (is.null(groups)) {
      stop("Paramter 'groups' cannot be NULL if type='grouping'!")
    }

    other_var <- unique(unlist(x$param$group))
    if (any(!groups %in% other_var)) {
      stop("Variables in 'groups' must contain variables from x$params$group!")
    }

    match_group <- lapply(
      x$param$group,
      function(z) {
        if (!is.null(z)) {
          identical(sort(z), sort(groups))
        }
      }
    )
    match_group <- any(unlist(match_group))
    if (length(groups) > 2) {
      stop("More than a maximum of 2 groups is not supported as of right now.")
    }
    if (!match_group) {
      if (length(groups) == 1) {
        stop("No results for ", groups, " present in the data")
      } else {
        stop("No results for the combination of ", groups[1], " and ",
             groups[2], " present in the data")
      }
    }
  }
  if (!sd.type %in% c("ribbon", "dot")) {
    stop("Parameter 'sd.type' can only take values 'ribbon' or 'dot'!")
  }
  # convert all factors for data tables in x
  # otherwise ther will be issues with merges
  x$Estimates <- convert_factors(x$Estimates)
  x$smallGroups <- convert_factors(x$smallGroups)
  x$cvHigh <- convert_factors(x$cvHigh)
  x$stEDecrease <- convert_factors(x$stEDecrease)
  #################
  # get variables from x and merge tables together
  period <- x$param$period

  poss.values <- c("Missing", "cvHigh", "SmallGroup+cvHigh", "SmallGroup", "OK")
  poss.color <- c("grey", "orangered3", "yellow2", "dodgerblue1", "green4")

  plot.x <- x$Estimates
  plot.x <- define_type(plot.x, x, variable = variable)
  values.plot <- plot.x[, unique(res_type)]
  color.plot <- poss.color[poss.values %in% values.plot]

  # define groups
  dt.eval("plot.x[grepl('^[[:digit:]]+$',", period,
          "),GROUPING:='Single Periods']")
  dt.eval("plot.x[grepl('-',", period, "),GROUPING:='Other']")
  dt.eval("plot.x[grepl('_',", period, "),GROUPING:='", x$param$period.mean,
          "-Period-Mean']")
  plot.x[, GROUPING := factor(GROUPING, levels = c("Single Periods", paste0(
    x$param$period.mean, "-Period-Mean"), "Other"))]

  #################
  # plot type 'summary'
  if (type == "summary") {
    title <- paste("Results of Groups per period for variable", variable)
    p1 <- ggplot(plot.x, aes(get(period), fill = res_type)) +
      geom_bar() + xlab("") + ylab("Count") + coord_flip() +
      facet_grid(GROUPING ~ ., scales = "free_y", space = "free_y") +
      theme(legend.title = element_blank()) +
      scale_fill_manual(breaks = values.plot,
                        values = color.plot) +
      ggtitle(title)
    plot(p1)
  } else if (type == "grouping") {
    #################
    # plot type 'grouping'
    if (length(groups) == 1) {
      groups <- c(period, groups)
    }

    if (length(other_var) > 0) {
      na_var <- other_var[!other_var %in% groups]
      if (length(na_var) > 0) {
        na_var <- paste(paste0("is.na(", na_var, ")"), collapse = "&")
        na_var <- paste0("&", na_var)
      }

      exp1 <- paste0("c(!is.na(", groups[1], ")&!is.na(", groups[2], ")",
                     na_var, ")")

      na_var <- other_var[!other_var %in% groups[1]]
      exp2_2 <- paste(paste0("is.na(", na_var, ")"), collapse = "&")
      exp2 <- paste0("c(!is.na(", groups[1], ")&", exp2_2, ")")
    } else {
      exp1 <- paste0("c(!is.na(", groups[1], ")&!is.na(", groups[2], "))")
      exp2 <- paste0("c(!is.na(", groups[1], ")&is.na(", groups[2], "))")
    }


    plot.group1 <- dt.eval("plot.x[", exp1, "]")
    plot.group2 <- dt.eval("plot.x[", exp2, "]")

    val_var <- paste0("val_", variable)
    ste_var <- paste0("stE_", variable)
    ste_bool <- variable


    # prepare plot.group1 for plotting
    plot.group1.period <- plot.group1[GROUPING == "Single Periods"]
    plot.group1.periodmean <- plot.group1[GROUPING == paste0(
      x$param$period.mean, "-Period-Mean")]

    # prepare plot.group2 for plotting
    plot.group2.period <- plot.group2[GROUPING == "Single Periods"]
    plot.group2.periodmean <- plot.group2[GROUPING == paste0(
      x$param$period.mean, "-Period-Mean")]

    # rename data columns for plot.group2
    ste_var2 <- paste0(ste_var, "2")
    val_var2 <- paste0(val_var, "2")
    ste_bool2 <- paste0(ste_bool, "2")

    # prepare output differently if period.mean has been used
    if (!is.null(x$param$period.mean)) {
      period.mean <- TRUE

      ste_var_mean <- paste0(ste_var, "_mean")
      ste_bool_mean <- paste0(ste_bool, "_mean")
      ste_var_mean2 <- paste0(ste_var_mean, "2")
      ste_bool_mean2 <- paste0(ste_bool_mean, "2")

      dt.eval("plot.group1.periodmean[,", period, ":=tstrsplit(", period,
              ",split='_',keep=", round((x$param$period.mean + 1) / 2), ")]")

      setnames(plot.group1.periodmean, c(ste_var, ste_bool),
               c(ste_var_mean, ste_bool_mean))

      plot.group1.period <- merge(
        plot.group1.period, dt.eval(
          "plot.group1.periodmean[,.(",
          paste(c(period, other_var, ste_var_mean, ste_bool_mean),
                collapse = ","),
          ")]"),
        all.x = TRUE)

      dt.eval("plot.group2.periodmean[,", period, ":=tstrsplit(", period,
              ",split='_',keep=", round((x$param$period.mean + 1) / 2), ")]")

      setnames(plot.group2.periodmean, c(ste_var, ste_bool),
               c(ste_var_mean, ste_bool_mean))
      plot.group2.period <- merge(
        plot.group2.period,
        dt.eval("plot.group2.periodmean[,.(",
                paste(c(period, other_var, ste_var_mean, ste_bool_mean),
                      collapse = ","), ")]"),
        all.x = TRUE)

      change.names <- c(ste_var, val_var, ste_bool, ste_var_mean, ste_bool_mean)
    } else {
      period.mean <- FALSE
      change.names <- c(ste_var, val_var, ste_bool)
    }

    setnames(plot.group2.period, change.names, paste0(change.names, "2"))
    select.group2 <- paste0(change.names, "2")

    # merge with plot.group
    if (period == groups[1]) {
      plot.group1.period <- merge(
        plot.group1.period,
        plot.group2.period[, mget(c(groups[1], select.group2))],
        all.x = TRUE, by = groups[1])
    } else {
      plot.group1.period <- merge(
        plot.group1.period,
        plot.group2.period[, mget(c(groups[1], period, select.group2))],
        all.x = TRUE, by = c(groups[1], period))
    }

    dt.eval("plot.group1.period[,", groups[1], ":=factor(", groups[1], ")]")
    dt.eval("plot.group1.period[,", groups[2], ":=factor(", groups[2], ")]")

    # plot results of subgroups and add results for mean over k periods
    p1 <- dt.eval("ggplot(plot.group1.period,aes(", period, ",", val_var,
                  ",group=", groups[2], "))+
                    geom_line(aes(colour=get(groups[2])))")

    if (sd.type == "ribbon") {
      p1 <- p1 + geom_ribbon(
        aes(ymin = get(val_var) - get(ste_var),
            ymax = get(val_var) + get(ste_var),
            fill = get(groups[2]), colour = get(groups[2])),
        linetype = 2, alpha = 0.1)


      # add results for group[1] -> higher level grouping
      p1 <- dt.eval("p1 + geom_line(aes(", period, ",", val_var2,
                    ",linetype='dotted'),color='black')")
      p1 <- p1 + geom_ribbon(
        aes(ymin = get(val_var2) - get(ste_var2),
            ymax = get(val_var2) + get(ste_var2)),
        fill = "grey", linetype = 2,  alpha = 0.1)

      if (period.mean) {
        # add results for k-mean-periods
        p1 <- p1 + geom_ribbon(
          aes(ymin = get(val_var) - get(ste_var_mean),
              ymax = get(val_var) + get(ste_var_mean), fill = get(groups[2]),
              colour = get(groups[2])),
          linetype = 2, alpha = 0.5)

        # add results for group[1] and k-period-mean
        p1 <- p1 + geom_ribbon(
          aes(ymin = get(val_var2) - get(ste_var_mean2),
              ymax = get(val_var2) + get(ste_var_mean2)),
          fill = "grey", linetype = 2, alpha = 0.5)
      }

    } else {
      plot.group1.period[, shape_bool := "Low st.Error"]
      plot.group1.period[get(ste_bool) == TRUE, shape_bool := "high st.Error"]
      plot.group1.period[, shape_bool2 := "Low st.Error"]
      plot.group1.period[get(ste_bool2) == TRUE, shape_bool2 := "high st.Error"]

      if (period.mean) {
        plot.group1.period[get(ste_bool_mean) == TRUE,
                           shape_bool := "high st.Error for mean"]
        plot.group1.period[get(ste_bool_mean2) == TRUE,
                           shape_bool2 := "high st.Error for mean"]
      }

      p1 <- dt.eval(
        "p1 + geom_point(data=plot.group1.period[get(ste_bool)==TRUE],
                         aes(", period, ",", val_var, ",shape=shape_bool))")

      # add results for group[1] -> higher level grouping
      p1 <- dt.eval("p1 + geom_line(aes(", period, ",", val_var2,
                    ",linetype='solid'))")
      p1 <- dt.eval(
        "p1 + geom_point(data=plot.group1.period[get(ste_bool2)==TRUE],
                         aes(", period, ",", val_var, ",shape=shape_bool2))")

      if (period.mean) {
        # add results for k-mean-periods
        p1 <- dt.eval("p1 + geom_point(data=plot.group1.period[", ste_bool_mean,
                      "==TRUE],
                      aes(", period, ",", val_var, ",shape=shape_bool))")

        # add results for group[1] and k-period-mean
        p1 <- dt.eval("p1 + geom_point(data=plot.group1.period[",
                      ste_bool_mean2, "==TRUE],
                      aes(", period, ",", val_var, ",shape=shape_bool2))")

      }
    }


    # define paramters for plot
    ylabel <- paste0("Point Estimate on variable ", x$param$var)
    p1 <- p1 + xlab(x$param$period) + ylab(ylabel) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))


    # split the plots into facets
    if (groups[1] == period) {
      p1 <- p1 + facet_wrap(~get(groups[2]))
      if (sd.type == "ribbon") {
        p1 <- p1 + guides(fill = 'none',
                          colour = 'none',
                          linetype = guide_legend(title = ""))
        p1 <- p1 + scale_linetype_manual(
          values = c("solid"), labels = paste0("Result for ", groups[1]))
      } else {
        shapes <- c(1, 16)
        names(shapes) <- c("high st.Error", "high st.Error for mean")
        p1 <- p1 + scale_shape_manual(values = shapes)
        p1 <- p1 + scale_colour_discrete(guide = 'none') +
        # p1 <- p1 + scale_colour_discrete(value="grey")
           theme(legend.title = element_blank())
        p1 <- p1 + scale_linetype_manual(values = c("dotted"), labels = paste0(
          "Result for ", groups[1]))
      }
    } else {
      p1 <- p1 + facet_wrap(~get(groups[1]))
      if (sd.type == "ribbon") {
        p1 <- p1 + guides(fill = guide_legend(title = groups[2]),
                          colour = guide_legend(title = groups[2]),
                          linetype = guide_legend(title = ""))
        p1 <- p1 + scale_linetype_manual(values = c("solid"), labels = paste0(
          "Result for ", groups[1]))
      } else {
        shapes <- c(1, 16)
        names(shapes) <- c("high st.Error", "high st.Error for mean")
        p1 <- p1 + scale_shape_manual(values = shapes) +
          guides(colour = guide_legend(title = groups[2]),
                 linetype = guide_legend(title = ""),
                 shape = guide_legend(title = ""))

        p1 <- p1 + scale_linetype_manual(values = c("dotted"), labels = paste0(
          "Result for ", groups[1]))
      }
    }

    plot(p1)
  }
}


define_type <- function(plot.x, x, variable = "HX080") {

  SMALLGROUP <- res_type <- NULL

  val_variable <- paste("val", variable, sep = "_")
  poss.values <- c("Missing", "cvHigh", "SmallGroup+cvHigh", "SmallGroup", "OK")

  # merge data
  merge.on <- copy(colnames(x$smallGroups))
  x$smallGroups[, SMALLGROUP := TRUE]
  plot.x <- x$smallGroups[plot.x,, on = c(merge.on)] # nolint
  merge.on <- unique(c(x$param$period, unlist(x$param$group)))
  plot.x <- x$cvHigh[plot.x,, on = c(merge.on)] # nolint

  # define type
  plot.x[!is.na(SMALLGROUP) & get(variable) == FALSE, res_type := "SmallGroup"]
  plot.x[!is.na(SMALLGROUP) & get(variable) == TRUE,
         res_type := "SmallGroup+cvHigh"]
  plot.x[is.na(get(val_variable)), res_type := "Missing"]
  plot.x[get(variable) == TRUE & is.na(res_type), res_type := "cvHigh"]
  plot.x[is.na(res_type), res_type := "OK"]
  plot.x[, res_type := factor(res_type, levels = poss.values)]

  return(plot.x)
}

convert_factors <- function(x) {
  x <- x[, lapply(
    .SD,
    function(z) {
      if (is.factor(z)) {
        as.character(z)
      } else {
        z
      }
    }
  )]
  return(x)
}
