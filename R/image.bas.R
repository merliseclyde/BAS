#' Images of models used in Bayesian model averaging
#'
#' Creates an image of the models selected using \code{\link{bas}}.
#'
#' Creates an image of the model space sampled using \code{\link{bas}}.  If a
#' subset of the top models are plotted, then probabilities are renormalized
#' over the subset.
#'
#'
#' @aliases image.bas image
#' @param x A BMA object of type 'bas' created by BAS
#' @param top.models Number of the top ranked models to plot
#' @param intensity Logical variable, when TRUE image intensity is proportional
#' to the probability or log(probability) of the model, when FALSE, intensity
#' is binary indicating just presence (light) or absence (dark) of a variable.
#' @param prob Logical variable for whether the area in the image for each
#' model should be proportional to the posterior probability (or log
#' probability) of the model (TRUE) or with equal area (FALSE).
#' @param log Logical variable indicating whether the intensities should be
#' based on log posterior odds (TRUE) or posterior probabilities (FALSE).  The
#' log of the posterior odds is for comparing the each model to the worst model
#' in the top.models.
#' @param rotate Should the image of models be rotated so that models are on
#' the y-axis and variables are on the x-axis (TRUE)
#' @param color The color scheme for image intensities. The value "rainbow"
#' uses the rainbow palette. The value "blackandwhite" produces a black and
#' white image (greyscale image)
#' @param subset indices of variables to include/exclude in plot
#' @param drop.always.included logical variable to drop variables that are
#' always forced into the model.  FALSE by default.
#' @param offset numeric value to add to intensity
#' @param digits number of digits in posterior probabilities to keep
#' @param vlas las parameter for placing variable names; see par
#' @param plas las parameter for posterior probability axis
#' @param rlas las parameter for model ranks
#' @param ... Other parameters to be passed to the \code{image} and \code{axis}
#' functions.
#' @note Suggestion to allow area of models be proportional to posterior
#' probability due to Thomas Lumley
#' @author Merlise Clyde \email{clyde@@stat.duke.edu}
#' @seealso \code{\link{bas}}
#' @references Clyde, M. (1999) Bayesian Model Averaging and Model Search
#' Strategies (with discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P.
#' Dawid, J.O. Berger, and A.F.M. Smith eds. Oxford University Press, pages
#' 157-185.
#' @keywords regression
#' @examples
#'
#' require(graphics)
#' data("Hald")
#' hald.ZSprior =  bas.lm(Y~ ., data=Hald,  prior="ZS-null")
#' image(hald.ZSprior, drop.always.included=TRUE) #drop the intercept
#'
#' @rdname image.bas
#' @family bas methods
#' @family bas plots
#' @method image bas
#' @export
image.bas <- function (x, top.models=20, intensity=TRUE, prob=TRUE, log=TRUE, rotate=TRUE, color="rainbow", subset=NULL, drop.always.included = FALSE,
                       offset=.75, digits=3, vlas=2,plas=0,rlas=0, ...)
{
  postprob = x$postprobs
  top.models = min(top.models, x$n.models)
  best = order(-x$postprobs)[1:top.models]
  postprob=postprob[best]/sum(postprob[best])
  which.mat <-  list2matrix.which(x, best)
  nvar <- ncol(which.mat)


  if (is.null(subset)) subset=1:nvar
  if (drop.always.included) {
    keep = x$include.always
    if (is.null(keep)) keep = 1
    subset = subset[!subset %in% keep]
    if (length(subset) == 0) stop("no models in subset to show; modify subset or drop.always.included")
  }

  which.mat =  which.mat[,subset, drop=FALSE]
  nvar = ncol(which.mat)
  namesx = x$namesx[subset]

  scale = postprob
  prob.lab= "Posterior Probability"

  if (log)   {
    scale = log(postprob) - min(log(postprob))
    prob.lab="Log Posterior Odds"
    # fix problem when scale has duplicate zeros
    zeros = which(scale == 0.0)
    nzeros = length(zeros)

    if (nzeros > 1) {
      scale[zeros] = seq(scale[zeros[1]-1],0.0,length=nzeros)/1000
    }
  }

  if (intensity)  which.mat = sweep(which.mat, 1, scale+offset,"*")

  if (rotate) scale = rev(scale)

  if (prob) m.scale = cumsum(c(0,scale))
  else  m.scale = seq(0, top.models)

  mat = (m.scale[-1] +  m.scale[-(top.models+1)])/2

  colors = switch(color,
    "rainbow" = c("black", rainbow(top.models+1, start=.75, end=.05)),
    "blackandwhite" =  gray(seq(0, 1, length=top.models))
    )

  par.old = par()$mar


  if (rotate) {
      par(mar = c(6,6,3,5) + .1)
      image(0:nvar, mat, t(which.mat[top.models:1, , drop=FALSE]),
            xaxt="n", yaxt="n",
            ylab="",
            xlab="",
            zlim=c(0, max(which.mat)),
            col=colors, ...)

      axis(2,at=mat,labels=round(scale, digits=digits), las=plas, ...)
      axis(4,at=mat,labels=top.models:1, las=rlas, ...)
      mtext("Model Rank", side=4, line=3, las=0)
      mtext(prob.lab, side=2, line=4, las=0)
      axis(1,at=(1:nvar -.5), labels=namesx, las=vlas, ...)
    }
  else{
    par(mar = c(6,8,6,2) + .1)
    image(mat, 0:nvar, which.mat[ , nvar:1,  drop=FALSE],
          xaxt="n", yaxt="n",
          xlab="",
          ylab="",
          zlim=c(0, max(which.mat)),
          col=colors, ...)

    axis(1,at=mat,labels=round(scale, digits=digits), las=plas,...)
    axis(3,at=mat,labels=1:top.models, las=rlas, ...)
    mtext("Model Rank", side=3, line=3)
    mtext(prob.lab, side=1, line=4)
    axis(2,at=(1:nvar -.5), labels=rev(namesx), las=vlas, ...)
  }

  box()

  par(par.old)
  invisible()
}
