# Protein Activity Data

This data sets includes several predictors of protein activity from an
experiment run at Glaxo.

## Format

`protein` is a dataframe with 96 observations and 8 predictor variables
of protein activity:

|         |           |         |                       |
|---------|-----------|---------|-----------------------|
| \[,1\]  | buf       | factor  | Buffer                |
| \[,2\]  | pH        | numeric |                       |
| \[,3\]  | NaCl      | numeric |                       |
| \[,4\]  | con       | numeric | protein concentration |
| \[,5\]  | ra        | factor  | reducing agent        |
| \[,6\]  | det       | factor  | detergent             |
| \[,7\]  | MgCl2     | numeric |                       |
| \[,8\]  | temp      | numeric | (temperature)         |
| \[,9\]  | prot.act1 | numeric |                       |
| \[,10\] | prot.act2 | numeric |                       |
| \[,11\] | prot.act3 | numeric |                       |
| \[,12\] | prot.act4 | numeric | protein activity      |

## Source

Clyde, M. A. and Parmigiani, G. (1998), Protein Construct Storage:
Bayesian Variable Selection and Prediction with Mixtures, Journal of
Biopharmaceutical Statistics, 8, 431-443
