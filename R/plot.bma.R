plot.bma = function (x, which = c(1:4),
  caption = c("Residuals vs Fitted", "Model Probabilities",
              "Model Complexity", "Inclusion Probabilities"),
  panel = if (add.smooth) panel.smooth else points, 
  sub.caption = NULL, main = "",
  ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...,
  id.n = 3,
  labels.id = names(residuals(x)), 
  cex.id = 0.75,  add.smooth = getOption("add.smooth"), 
  label.pos = c(4, 2))
{
    if (!inherits(x, "bma")) 
      stop("use only with \"bma\" objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
      stop("'which' must be in 1:3")
    show <- rep(FALSE, 2)
    show[which] <- TRUE

    yhat = fitted(x, type="BMA")
    r = x$Y - yhat
    n <- length(r)
    if (id.n > 0) {
      if (is.null(labels.id)) 
        labels.id <- paste(1:n)
      iid <- 1:id.n
      show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
      text.id <- function(x, y, ind, adj.x = TRUE) {
        labpos <- if (adj.x) 
          label.pos[1 + as.numeric(x > mean(range(x)))]
        else 3
        text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE, 
             pos = labpos, offset = 0.25)
      }
      
      if (any(show[2:3])) {
        show.m = sort.list(x$logmarg, decreasing = TRUE)[iid]
        label.m = paste(1:x$n.models)
      }
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
      else cc[1]
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
      plot(yhat, r, xlab = "Predictions under BMA",
           ylab = "Residuals", main = main, 
           ylim = ylim, type = "n", ...)
      panel(yhat, r, ...)
      if (one.fig) 
        title(sub = sub.caption, ...)
      mtext(caption[1], 3, 0.25)
      if (id.n > 0) {
        y.id <- r[show.r]
        y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
        text.id(yhat[show.r], y.id, show.r)
      }
      abline(h = 0, lty = 3, col = "gray")
    }
    if (show[2]) {
      cum.prob = cumsum(x$postprobs)
      m.index = 1:x$n.models
      ylim <- range(cum.prob, na.rm = TRUE)
      ylim[2] <- ylim[2] + diff(ylim) * 0.075
      plot(m.index, cum.prob,
           xlab="Model Search Order", ylab="Cumulative Probability",
           type="n")
      panel(m.index,cum.prob)
      if (one.fig) 
        title(sub = sub.caption, ...)
      mtext(caption[2], 3, 0.25)
      if (id.n > 0) 
        text.id(m.index[show.m], cum.prob[show.m], show.m)
    }
    if (show[3]) {
      logmarg = x$logmarg
      dim = x$size
      ylim <- range(logmarg, na.rm = TRUE)
      plot(dim, logmarg,
           xlab = "Model Dimension", ylab = "log(Marginal)",
           main = main, 
           ylim = ylim, ...)
      if (one.fig) 
        title(sub = sub.caption, ...)
      mtext(caption[3], 3, 0.25)
      if (id.n > 0) 
        text.id(dim[show.m], logmarg[show.m], show.m)
    }
    if (show[4]) {
      probne0 = x$probne0
      variables = 1:x$n.vars
      ylim <- c(0,1)
      plot(variables, probne0,
           xlab = "", ylab = "Marginal Inclusion Probability",
           xaxt="n",
           main = main, type="h", col=(1+(probne0>= .5)),
           ylim = ylim, ...)
      if (one.fig) 
        title(sub = sub.caption, ...)
      mtext(x$namesx, side=1, line=0.25, at=variables, las=2, cex=.75 )
      mtext(caption[4], 3, 0.25)
      if (id.n > 0) 
        text.id(dim[show.m], logmarg[show.m], show.m)
    }
    
    if (!one.fig && par("oma")[3] >= 1) {
      mtext(sub.caption, outer = TRUE, cex = 1.25)}
    invisible()
    
  }  
