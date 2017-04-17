

#' Bayesian Model Averaging using Bayesian Adaptive Sampling
#' 
#' Package for Bayesian Model Averaging in linear models using stochastic or
#' deterministic sampling without replacement from posterior distributions.
#' Prior distributions on coefficients are of the form of Zellner's g-prior or
#' mixtures of g-priors. Options include the Zellner-Siow Cauchy Priors, the
#' Liang et al hyper-g priors, Local and Global Empirical Bayes estimates of g,
#' and other default model selection criteria such as AIC and BIC. Sampling
#' probabilities may be updated based on the sampled models.
#' 
#' \tabular{ll}{ Package: \tab BAS\cr Depends: \tab R (>= 2.8)\cr License: \tab
#' GPL-2\cr URL: http://www.stat.duke.edu/~clyde\cr }
#' 
#' Index: \preformatted{ }
#' 
#' @name BAS-package
#' @aliases BAS-package BAS
#' @docType package
#' @author Merlise Clyde, \cr Maintainer: Merlise Clyde <clyde@@stat.duke.edu>
#' @seealso \code{\link[BAS]{bas}}
#' @references Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive
#' Sampling for Variable Selection and Model Averaging. Journal of
#' Computational Graphics and Statistics.  20:80-101 \cr
#' \url{http://dx.doi.org/10.1198/jcgs.2010.09049}
#' 
#' Clyde, M. and George, E. I. (2004) Model uncertainty. Statist. Sci., 19,
#' 81-94. \cr \url{http://dx.doi.org/10.1214/088342304000000035}
#' 
#' Clyde, M. (1999) Bayesian Model Averaging and Model Search Strategies (with
#' discussion). In Bayesian Statistics 6. J.M. Bernardo, A.P. Dawid, J.O.
#' Berger, and A.F.M. Smith eds. Oxford University Press, pages 157-185.
#' 
#' Li, Y. and Clyde, M. (2015) Mixtures of g-priors in Generalized Linear
#' Models.  \url{http://arxiv.org/abs/1503.06913}
#' 
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J.O. (2008) Mixtures
#' of g-priors for Bayesian Variable Selection. Journal of the American
#' Statistical Association. 103:410-423.  \cr
#' \url{http://dx.doi.org/10.1198/016214507000001337}
#' @keywords package regression
#' @examples
#' 
#' demo(BAS.USCrime)
#' demo(BAS.hald)
#' 
NULL





#' Bodyfat Data
#' 
#' Lists estimates of the percentage of body fat determined by underwater
#' weighing and various body circumference measurements for 252 men.  Accurate
#' measuremnt of body fat is inconvenient/costly and it is desirable to have
#' easy methods of estimating body fat that are not inconvenient/costly.
#' 
#' A variety of popular health books suggest that the readers assess their
#' health, at least in part, by estimating their percentage of body fat. In
#' Bailey (1994), for instance, the reader can estimate body fat from tables
#' using their age and various skin-fold measurements obtained by using a
#' caliper. Other texts give predictive equations for body fat using body
#' circumference measurements (e.g. abdominal circumference) and/or skin-fold
#' measurements. See, for instance, Behnke and Wilmore (1974), pp. 66-67;
#' Wilmore (1976), p. 247; or Katch and McArdle (1977), pp. 120-132).#
#' 
#' Percentage of body fat for an individual can be estimated once body density
#' has been determined. Folks (e.g. Siri (1956)) assume that the body consists
#' of two components - lean body tissue and fat tissue. Letting
#' 
#' D = Body Density (gm/cm^3) A = proportion of lean body tissue B = proportion
#' of fat tissue (A+B=1) a = density of lean body tissue (gm/cm^3) b = density
#' of fat tissue (gm/cm^3)
#' 
#' we have D = 1/[(A/a) + (B/b)] and solving for B we find B = (1/D)*[ab/(a-b)]
#' - [b/(a-b)].
#' 
#' Using the estimates a=1.10 gm/cm^3 and b=0.90 gm/cm^3 (see Katch and McArdle
#' (1977), p. 111 or Wilmore (1976), p. 123) we come up with "Siri's equation":
#' 
#' Percentage of Body Fat (i.e. 100*B) = 495/D - 450.#
#' 
#' Volume, and hence body density, can be accurately measured a variety of
#' ways. The technique of underwater weighing "computes body volume as the
#' difference between body weight measured in air and weight measured during
#' water submersion. In other words, body volume is equal to the loss of weight
#' in water with the appropriate temperature correction for the water's
#' density" (Katch and McArdle (1977), p. 113). Using this technique,
#' 
#' Body Density = WA/[(WA-WW)/c.f. - LV]
#' 
#' where WA = Weight in air (kg) WW = Weight in water (kg) c.f. = Water
#' correction factor (=1 at 39.2 deg F as one-gram of water occupies exactly
#' one cm^3 at this temperature, =.997 at 76-78 deg F) LV = Residual Lung
#' Volume (liters)
#' 
#' (Katch and McArdle (1977), p. 115). Other methods of determining body volume
#' are given in Behnke and Wilmore (1974), p. 22 ff.
#' 
#' Measurement standards are apparently those listed in Benhke and Wilmore
#' (1974), pp. 45-48 where, for instance, the abdomen circumference is measured
#' "laterally, at the level of the iliac crests, and anteriorly, at the
#' umbilicus".)
#' 
#' @name bodyfat
#' @aliases Bodyfat bodyfat
#' @docType data
#' @format A data frame with 252 observations on the following 15 variables.
#' \describe{ \item{list("Density")}{a numeric vector for the density
#' determined from underwater weighing} \item{list("Bodyfat")}{percent body fat
#' from Siri's (1956) equation} \item{list("Age")}{age of individual in years}
#' \item{list("Weight")}{weight of the individual in pounds}
#' \item{list("Height")}{height of individual in inches}
#' \item{list("Neck")}{neck circumference in centimeters (cm)}
#' \item{list("Chest")}{chest circumference (cm)}
#' \item{list("Abdomen")}{abdomen circumference (cm)} \item{list("Hip")}{hip
#' circumference (cm)} \item{list("Thigh")}{thigh circumference (cm)}
#' \item{list("Knee")}{knee circumference (cm)} \item{list("Ankle")}{ankle
#' circumference (cm)} \item{list("Biceps")}{bicep (extended) circumference
#' (cm)} \item{list("Forearm")}{forearm circumference (cm)}
#' \item{list("Wrist")}{wrist circumference (cm)} }
#' @references Bailey, Covert (1994). Smart Exercise: Burning Fat, Getting Fit,
#' Houghton-Mifflin Co., Boston, pp. 179-186.
#' 
#' Behnke, A.R. and Wilmore, J.H. (1974). Evaluation and Regulation of Body
#' Build and Composition, Prentice-Hall, Englewood Cliffs, N.J.
#' 
#' Siri, W.E. (1956), "Gross composition of the body", in Advances in
#' Biological and Medical Physics, vol. IV, edited by J.H. Lawrence and C.A.
#' Tobias, Academic Press, Inc., New York.
#' 
#' Katch, Frank and McArdle, William (1977). Nutrition, Weight Control, and
#' Exercise, Houghton Mifflin Co., Boston.
#' 
#' Wilmore, Jack (1976). Athletic Training and Physical Fitness: Physiological
#' Principles of the Conditioning Process, Allyn and Bacon, Inc., Boston.
#' @source These data are used to produce the predictive equations for lean
#' body weight given in the abstract "Generalized body composition prediction
#' equation for men using simple measurement techniques", K.W. Penrose, A.G.
#' Nelson, A.G. Fisher, FACSM, Human Performance Research Center, Brigham Young
#' University, Provo, Utah 84602 as listed in _Medicine and Science in Sports
#' and Exercise_, vol. 17, no. 2, April 1985, p. 189. (The predictive equations
#' were obtained from the first 143 of the 252 cases that are listed below).
#' The data were generously supplied by Dr. A. Garth Fisher who gave permission
#' to freely distribute the data and use for non-commercial purposes.
#' @keywords datasets
#' @examples
#' 
#' data(bodyfat)
#' bodyfat.bas = bas.lm(Bodyfat ~ Abdomen, data=bodyfat, prior="ZS-null")
#' summary(bodyfat.bas)
#' plot(Bodyfat ~ Abdomen, data=bodyfat, xlab="abdomen circumference (cm)")
#' betas = coef(bodyfat.bas)$postmean   # current version has that intercept is ybar
#' betas[1] = betas[1] - betas[2]*bodyfat.bas$mean.x
#' abline(betas)
#' abline(coef(lm(Bodyfat ~ Abdomen, data=bodyfat)), col=2, lty=2)
#' 
NULL









#' Hald Data
#' 
#' The Hald data have been used in many books and papers to illustrate variable
#' selection. The data relate to an engineering application that was concerned
#' with the effect of the composition of cement on heat evolved during
#' hardening. The response variable \emph{Y} is the \emph{heat evolved} in a
#' cement mix. The four explanatory variables are ingredients of the mix, X1:
#' \emph{tricalcium aluminate}, X2: \emph{tricalcium silicate}, X3:
#' \emph{tetracalcium alumino ferrite}, X4: \emph{dicalcium silicate}. An
#' important feature of these data is that the variables X1 and X3 are highly
#' correlated, as well as the variables X2 and X4.  Thus we should expect any
#' subset of (X1,X2,X3,X4) that includes one variable from highly correlated
#' pair to do as any subset that also includes the other member.
#' 
#' 
#' @name Hald
#' @aliases Hald hald
#' @docType data
#' @format \code{hald} is a dataframe with 13 observations and 5 variables
#' (columns),
#' 
#' Y: Heat evolved per gram of cement (in calories) X1: Amount of tricalcium
#' aluminate X2: Amount of tricalcium silicate X3: Amount of tetracalcium
#' alumino ferrite X4: Amount of dicalcium silicate
#' @source Wood, H., Steinour, H.H., and Starke, H.R. (1932). "Effect of
#' Composition of Portland cement on Heat Evolved During Hardening", Industrila
#' and Engineering Chemistry, 24, 1207-1214.
#' @keywords datasets
NULL





#' Protein Activity Data
#' 
#' This data sets includes several predictors of protein activity from an
#' experiment run at Glaxo.
#' 
#' 
#' @name protein
#' @docType data
#' @format \code{protein} is a dataframe with 96 observations and 8 predictor
#' variablesof protein activity: \tabular{llll}{ [,1] \tab buf \tab factor \tab
#' Buffer \cr [,2] \tab pH \tab numeric \tab \cr [,3] \tab NaCl \tab numeric
#' \tab \cr [,4] \tab con \tab numeric \tab protein concentration\cr [,5] \tab
#' ra \tab factor \tab reducing agent\cr [,6] \tab det \tab factor \tab
#' detergent\cr [,7] \tab MgCl2 \tab numeric\tab \cr [,8] \tab temp \tab
#' numeric\tab (temerature)\cr [,9] \tab prot.act1 \tab numeric\tab \cr [,10]
#' \tab prot.act2 \tab numeric\tab \cr [,11] \tab prot.act3 \tab numeric\tab
#' \cr [,12] \tab prot.act4 \tab numeric\tab protein activity }
#' @source Clyde, M. A. and Parmigiani, G. (1998), Protein Construct Storage:
#' Bayesian Variable Selection and Prediction with Mixtures, Journal of
#' Biopharmaceutical Statistics, 8, 431-443
#' @keywords datasets
NULL



