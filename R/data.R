#' Bodyfat Data
#'
#' Lists estimates of the percentage of body fat determined by underwater
#' weighing and various body circumference measurements for 252 men.  Accurate
#' measurement of body fat is inconvenient/costly and it is desirable to have
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
#' Measurement standards are apparently those listed in Behnke and Wilmore
#' (1974), pp. 45-48 where, for instance, the abdomen circumference is measured
#' "laterally, at the level of the iliac crests, and anteriorly, at the
#' umbilicus".)
#'
#' @name bodyfat
#' @aliases Bodyfat bodyfat
#' @docType data
#' @format A data frame with 252 observations on the following 15 variables.
#' \describe{ \item{Density}{a numeric vector for the density
#' determined from underwater weighing} \item{Bodyfat}{percent body fat
#' from Siri's (1956) equation} \item{Age}{age of individual in years}
#' \item{Weight}{weight of the individual in pounds}
#' \item{Height}{height of individual in inches}
#' \item{Neck}{neck circumference in centimeters (cm)}
#' \item{Chest}{chest circumference (cm)}
#' \item{Abdomen}{abdomen circumference (cm)} \item{Hip}{hip
#' circumference (cm)} \item{"Thigh"}{thigh circumference (cm)}
#' \item{"Knee"}{knee circumference (cm)} \item{Ankle}{ankle
#' circumference (cm)} \item{Biceps}{bicep (extended) circumference
#' (cm)} \item{Forearm}{forearm circumference (cm)}
#' \item{Wrist}{wrist circumference (cm)} }
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

#' Climate Data
#' @name climate
#' @docType data
#' @format Scientists are interested in the Earth's temperature change since the last
#' glacial maximum, about 20,000 years ago. The first study to estimate the
#' temperature change was published in 1980, and estimated a change of -1.5 degrees
#'  C, +/- 1.2 degrees C in tropical sea surface temperatures.
#'  The negative value means that the Earth was colder then than now.
#'  Since 1980 there have been many other studies.
#' \code{climate} is a dataset with 63 measurements on 5 variables:
#' \describe{\item{\emph{deltaT}}{ the response variables, which is the change in temperature
#' in degrees Celsius;}
#' \item{\emph{sdev}}{a standard deviation for the calculated \emph{deltaT};}
#' \item{\emph{proxy}}{a number 1-8 reflecting which type of measurement system was used to derive
#' deltaT. Some proxies can be used over land, others over water.
#' The proxies are coded as\cr
#' 1 "Mg/Ca"         \cr
#' 2 "alkenone"      \cr
#' 3 "Faunal"        \cr
#' 4 "Sr/Ca"         \cr
#' 5 "del 180"       \cr
#' 6 "Ice Core"      \cr
#' 7 "Pollen"        \cr
#' 8 "Noble Gas"     \cr
#'}
#'\item{\emph{T/M}}{, an indicator of whether it was a terrestrial or marine study (T/M),
#'  which is coded as 0 for Terrestrial, 1 for Marine;}
#'\item{ \emph{latitude}}{the latitude where the data were collected.}}
#' @source Data provided originally by Michael Lavine and available at \url{https://stat.duke.edu/sites/stat.duke.edu/files/climate.dat}
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
#' Composition of Portland cement on Heat Evolved During Hardening", Industrial
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
#' variables of protein activity: \tabular{llll}{ [,1] \tab buf \tab factor \tab
#' Buffer \cr [,2] \tab pH \tab numeric \tab \cr [,3] \tab NaCl \tab numeric
#' \tab \cr [,4] \tab con \tab numeric \tab protein concentration\cr [,5] \tab
#' ra \tab factor \tab reducing agent\cr [,6] \tab det \tab factor \tab
#' detergent\cr [,7] \tab MgCl2 \tab numeric\tab \cr [,8] \tab temp \tab
#' numeric\tab (temperature)\cr [,9] \tab prot.act1 \tab numeric\tab \cr [,10]
#' \tab prot.act2 \tab numeric\tab \cr [,11] \tab prot.act3 \tab numeric\tab
#' \cr [,12] \tab prot.act4 \tab numeric\tab protein activity }
#' @source Clyde, M. A. and Parmigiani, G. (1998), Protein Construct Storage:
#' Bayesian Variable Selection and Prediction with Mixtures, Journal of
#' Biopharmaceutical Statistics, 8, 431-443
#' @keywords datasets
#'
NULL

#'  Horseshoe Crab Data
#'
#'  Data on horseshoe crabs (\emph{Limulus polyphemus}).
#'  Responseis number of males surrounding a breeding female,
#'  color (factor), condition (factor), weight (quantitative),
#'  and width (quantitative) of the female.
#'
#' @name crabs
#' @docType data
#'
#' @format A data frame with 173 observations on 6 variables.
#'Individuals (rows of the data frame) are female horseshoe crabs.
#'Variables other than \code{satell} refer to these females. The variables are
#'  \describe{
#'    \item{color}{color. The colors given in
#'      Agresti are \dQuote{light medium}, \dQuote{medium}, \dQuote{dark medium},
#'      and \dQuote{dark}.  Here they are abbreviated to \code{light},
#'      \code{medium}, \code{dark}, and \code{darker}, respectively.}
#'    \item{spine}{spine condition.  The conditions given in Agresti are
#'      \dQuote{both good}, \dQuote{one worn or broken}, and
#'      \dQuote{both worn or broken}.
#'      Here they are abbreviated to \code{good}, \code{middle}, \code{bad},
#'      respectively.}
#'    \item{width}{carapace width in centimeters}
#'    \item{satell}{number of satellites, which males clustering around the
#'      female in addition to the male with which she is breeding.}
#'    \item{weight}{weight in grams.}
#'    \item{y}{shorthand for \code{as.numeric(satell > 0)}.}
#'  }
#' @details
#'  Quoting from the abstract of Brockmann (1996). \dQuote{Horseshoe crabs
#'    arrive on the beach in pairs and spawn \ldots during \ldots high tides.
#'    Unattached males also come to the beach, crowd around the nesting couples
#'    and compete with attached males for fertilizations. Satellite males form
#'    large groups around some couples while ignoring others, resulting in
#'    a nonrandom distribution that cannot be explained by local environmental
#'    conditions or habitat selection.}
#' @source Agresti, A. (2013)  \emph{Categorical Data Analysis},
#' Wiley, Hoboken, NJ., Section 4.3.2,
#'  \url{http://www.stat.ufl.edu/~aa/cda/data.html}
#'
#'  Brockmann, H. J. (1996)
#'  Satellite Male Groups in Horseshoe Crabs, \emph{Limulus polyphemus},
#'  \emph{Ethology}, \bold{102}, 1--21.
#' @keywords datasets
NULL
