# Hald Data

The Hald data have been used in many books and papers to illustrate
variable selection. The data relate to an engineering application that
was concerned with the effect of the composition of cement on heat
evolved during hardening. The response variable *Y* is the *heat
evolved* in a cement mix. The four explanatory variables are ingredients
of the mix, X1: *tricalcium aluminate*, X2: *tricalcium silicate*, X3:
*tetracalcium alumino ferrite*, X4: *dicalcium silicate*. An important
feature of these data is that the variables X1 and X3 are highly
correlated, as well as the variables X2 and X4. Thus we should expect
any subset of (X1,X2,X3,X4) that includes one variable from highly
correlated pair to do as any subset that also includes the other member.

## Format

`hald` is a dataframe with 13 observations and 5 variables (columns),

Y: Heat evolved per gram of cement (in calories) X1: Amount of
tricalcium aluminate X2: Amount of tricalcium silicate X3: Amount of
tetracalcium alumino ferrite X4: Amount of dicalcium silicate

## Source

Wood, H., Steinour, H.H., and Starke, H.R. (1932). "Effect of
Composition of Portland cement on Heat Evolved During Hardening",
Industrial and Engineering Chemistry, 24, 1207-1214.
