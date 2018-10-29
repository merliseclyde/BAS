#' Plot Diagnostics for an BAS Object
#'
#' Four plots (selectable by 'which') are currently available: a plot of
#' residuals against fitted values, Cumulative Model Probabilities, log
#' marginal likelihoods versus model dimension, and marginal inclusion
#' probabilities.
#'
#' This provides a panel of 4 plots: the first is a plot of the residuals
#' versus fitted values under BMA. The second is a plot of the cumulative
#' marginal likelihoods of models; if the model space cannot be enumerated then
#' this provides some indication of whether the probabilities are leveling off.
#' The third is a plot of log marginal likelihood versus model dimension and
#' the fourth plot show the posterior marginal inclusion probabilities.
#'
#' @param x \code{bas} BMA object result of 'bas'
#' @param which if a subset of the plots is required, specify a subset of the
#' numbers '1:4'
#' @param caption captions to appear above the plots
#' @param panel panel function.  The useful alternative to 'points',
#' 'panel.smooth' can be chosen by 'add.smooth = TRUE'
#' @param sub.caption common title-above figures if there are multiple; used as
#' 'sub' (s.'title') otherwise.  If 'NULL', as by default, a possible shortened
#' version of \code{deparse(x$call)} is used
#' @param main title to each plot-in addition to the above 'caption'
#' @param ask logical; if 'TRUE', the user is asked before each plot, see
#' 'par(ask=.)'
#' @param col.in color for the included variables
#' @param col.ex color for the excluded variables
#' @param col.pch color for points in panels 1-3
#' @param cex.lab graphics parameter to control size of variable names
#' @param ...  other parameters to be passed through to plotting functions
#' @param id.n number of points to be labeled in each plot, starting with the
#' most extreme
#' @param labels.id vector of labels, from which the labels for extreme points
#' will be chosen.  'NULL' uses observation numbers
#' @param cex.id magnification of point labels.
#' @param add.smooth logical indicating if a smoother should be added to most
#' plots; see also 'panel' above
#' @param label.pos positioning of labels, for the left half and right half of
#' the graph respectively, for plots 1-4
#' @param subset indices of variables to include/exclude in plot of marginal posterior
#' inclusion probabilities (NULL).
#' @param drop.always.included logical variable to drop marginal posterior inclusion
#' probabilities
#' for variables that are always forced into the model.  FALSE by default.
#' @author Merlise Clyde, based on plot.lm by John Maindonald and Martin
#' Maechler
#' @seealso \code{\link{plot.coef.bas}} and \code{\link{image.bas}}.
#' @keywords regression
#' @examples
#'
#' data(Hald)
#' hald.gprior =  bas.lm(Y~ ., data=Hald, prior="g-prior", alpha=13,
#'                       modelprior=beta.binomial(1,1),
#'                       initprobs="eplogp")
#'
#' plot(hald.gprior)
#'
#'
#' @rdname plot
#' @family bas plots
#' @export
plot.bas = function (x,
                     which = c(1:4),
                     caption = c(
                       "Residuals vs Fitted",
                       "Model Probabilities",
                       "Model Complexity",
                       "Inclusion Probabilities"
                     ),
                     panel = if (add.smooth)
                       panel.smooth
                     else
                       points,
                     sub.caption = NULL,
                     main = "",
                     ask = prod(par("mfcol")) < length(which) && dev.interactive(),
                     col.in = 2,
                     col.ex = 1,
                     col.pch = 1,
                     cex.lab = 1,
                     ...,
                     id.n = 3,
                     labels.id = NULL,
                     cex.id = 0.75,
                     add.smooth = getOption("add.smooth"),
                     label.pos = c(4, 2),
                     subset = NULL,
                     drop.always.included = FALSE)
{
  if (!inherits(x, "bas"))
    stop("use only with \"bas\" objects")
  if (!is.numeric(which) || any(which < 1) || any(which > 4))
    stop("'which' must be in 1:4")
  show <- rep(FALSE, 4)
  show[which] <- TRUE

  iid <- 1:id.n

  if (show[1]) {
    yhat = fitted(x, estimator = "BMA")
    r = x$Y - yhat
    n <- length(r)
    if (id.n > 0) {
      if (is.null(labels.id))
        labels.id <- paste(1:n)
      show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
    }
  }

  text.id <- function(x, y, ind, adj.x = TRUE) {
    labpos <- if (adj.x)
      label.pos[1 + as.numeric(x > mean(range(x)))]
    else
      3
    text(
      x,
      y,
      labels.id[ind],
      cex = cex.id,
      xpd = TRUE,
      pos = labpos,
      offset = 0.25
    )
  }

  if (any(show[2:3])) {
    show.m = sort.list(x$logmarg, decreasing = TRUE)[iid]
    label.m = paste(1:x$n.models)
  }


  if (is.null(sub.caption)) {
    cal <- x$call
    if (!is.na(m.f <- match("formula", names(cal)))) {
      cal <- cal[c(1, m.f)]
      names(cal)[2] <- ""
    }
    cc <- deparse(cal, 80)
    nc <- nchar(cc[1])
    abbr <- length(cc) > 1 || nc > 75
    sub.caption <- if (abbr)
      paste(substr(cc[1], 1, min(75, nc)), "...")
    else
      cc[1]
  }

  one.fig <- prod(par("mfcol")) == 1
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }

  if (show[1]) {
    ylim <- range(r, na.rm = TRUE)
    if (id.n > 0)
      ylim <- extendrange(r = ylim, f = 0.08)
    plot(
      yhat,
      r,
      xlab = "Predictions under BMA",
      ylab = "Residuals",
      main = main,
      ylim = ylim,
      type = "n",
      col = col.pch,
      ...
    )
    panel(yhat, r, ...)
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(caption[1], 3, 0.25)
    if (id.n > 0) {
      y.id <- r[show.r]
      y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ") / 3
      text.id(yhat[show.r], y.id, show.r)
    }
    abline(h = 0, lty = 3, col = "gray")
  }
  if (show[2]) {
    cum.prob = cumsum(x$postprobs)
    m.index = 1:x$n.models
    ylim <- range(cum.prob, na.rm = TRUE)
    ylim[2] <- ylim[2] + diff(ylim) * 0.075
    plot(
      m.index,
      cum.prob,
      xlab = "Model Search Order",
      ylab = "Cumulative Probability",
      type = "n",
      col = col.pch,
      ...
    )
    panel(m.index, cum.prob)
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(caption[2], 3, 0.25)
    #if (id.n > 0)
    #  text.id(m.index[show.m], cum.prob[show.m], show.m)
  }
  if (show[3]) {
    logmarg = x$logmarg
    dim = x$size
    ylim <- range(logmarg, na.rm = TRUE)
    plot(
      dim,
      logmarg,
      xlab = "Model Dimension",
      ylab = "log(Marginal)",
      main = main,
      ylim = ylim,
      col = col.pch,
      ...
    )
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(caption[3], 3, 0.25)
    if (id.n > 0)
      text.id(dim[show.m], logmarg[show.m], show.m)
  }
  if (show[4]) {
    if (is.null(subset))
      subset = 1:x$n.vars
    if (drop.always.included) {
      keep = x$include.always
      if (is.null(keep))
        keep = 1
      subset = subset[!subset %in% keep]
      if (length(subset) == 0)
        stop("no models in subset to show; modify subset or drop.always.included")
    }
    probne0 = x$probne0[subset]
    nvars = length(subset)
    variables = 1:nvars
    ylim <- c(0, 1)
    colors = rep(0, nvars)
    colors[probne0 > .5] = col.in
    colors[probne0 <= .5] = col.ex

    plot(
      variables,
      probne0,
      xlab = "",
      ylab = "Marginal Inclusion Probability",
      xaxt = "n",
      main = main,
      type = "h",
      col = colors,
      ylim = ylim,
      ...
    )
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(
      x$namesx[subset],
      side = 1,
      line = 0.25,
      at = variables,
      las = 2,
      cex = cex.lab,
      ...
    )
    mtext(caption[4], 3, 0.25)
    #if (id.n > 0)
    # text.id(dim[show.m], logmarg[show.m], show.m)
  }

  if (!one.fig && par("oma")[3] >= 1) {
    mtext(sub.caption, outer = TRUE, cex = 1.25)
  }
  invisible()

}
