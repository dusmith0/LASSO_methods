# Coordinate-descent algorithm for LASSO

## Introduction

This was created to accomplish an assignment for my Reproducible Computation Course. The following is an overview of what I was expected to complete. 

We consider the training data consisting of
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
samples
![(x_i, y_i)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28x_i%2C%20y_i%29 "(x_i, y_i)"),
![x_i\in \mathbb{R}^p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_i%5Cin%20%5Cmathbb%7BR%7D%5Ep "x_i\in \mathbb{R}^p")
(vector of covariates for sample
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i")),
![y_i\in \mathbb{R}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_i%5Cin%20%5Cmathbb%7BR%7D "y_i\in \mathbb{R}")
(response for sample
![i](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;i "i"))
supplied as matrix
![X \in \mathbb{R}^{n \times p}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bn%20%5Ctimes%20p%7D "X \in \mathbb{R}^{n \times p}")
and vector
![Y\in\mathbb{R}^{n}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%5Cin%5Cmathbb%7BR%7D%5E%7Bn%7D "Y\in\mathbb{R}^{n}"),
respectively. We would like to fit a linear model

![Y = \beta_0 + X\beta + \varepsilon,](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%20%3D%20%5Cbeta_0%20%2B%20X%5Cbeta%20%2B%20%5Cvarepsilon%2C "Y = \beta_0 + X\beta + \varepsilon,")

where the sample size
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
is small compared to the number of covariates
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p").
I used LASSO algorithm to fit this model (find
![\beta_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_0 "\beta_0")
and
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")),
and use 5-fold cross-validation to select the tuning parameter.


## Part 1 - Implement LASSO with cross validation

For this part I implemented LASSO and 5 fold cross validation to select the tuning parameter in R.

1)  We first center
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y"),
    and center and scale
    ![X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X "X")
    to form
    ![\widetilde Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidetilde%20Y "\widetilde Y")
    and
    ![\widetilde X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidetilde%20X "\widetilde X"),
    and fit

    ![\widetilde Y = \widetilde X \widetilde \beta + \varepsilon.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidetilde%20Y%20%3D%20%5Cwidetilde%20X%20%5Cwidetilde%20%5Cbeta%20%2B%20%5Cvarepsilon. "\widetilde Y = \widetilde X \widetilde \beta + \varepsilon.")

    Compared to original model, there is no intercept
    ![\beta_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_0 "\beta_0")
    (because of centering), and
    ![\widetilde X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidetilde%20X "\widetilde X")
    is such that each column satisfies
    ![n^{-1}\widetilde X_j^{\top}\widetilde X_j = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%5E%7B-1%7D%5Cwidetilde%20X_j%5E%7B%5Ctop%7D%5Cwidetilde%20X_j%20%3D%201 "n^{-1}\widetilde X_j^{\top}\widetilde X_j = 1")
    (because of scaling).

2)  Next I solved the following LASSO problem for various
    ![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
    values

    ![\widetilde \beta = \arg\min\_{\beta}\left\\{(2n)^{-1}\\\|\widetilde Y-\widetilde X\beta\\\|\_2^2 + \lambda \\\|\beta\\\|\_1\right\\}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidetilde%20%5Cbeta%20%3D%20%5Carg%5Cmin_%7B%5Cbeta%7D%5Cleft%5C%7B%282n%29%5E%7B-1%7D%5C%7C%5Cwidetilde%20Y-%5Cwidetilde%20X%5Cbeta%5C%7C_2%5E2%20%2B%20%5Clambda%20%5C%7C%5Cbeta%5C%7C_1%5Cright%5C%7D. "\widetilde \beta = \arg\min_{\beta}\left\{(2n)^{-1}\|\widetilde Y-\widetilde X\beta\|_2^2 + \lambda \|\beta\|_1\right\}.")

    To solve LASSO, we will use coordinate-descent algorithm with **warm
    starts**.

3)  I used the
    ![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")-fold
    cross-validation to select
    ![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda"),
    and then find
    ![\beta_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_0 "\beta_0"),
    ![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
    for original
    ![X](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X "X")
    and
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    based on back-scaling and back-centering.
    

## Application to Riboflavin data

Your implementation will be used for analysis of riboflavin data
available from the R package . The **RiboflavinDataAnalysis.R** gives
starter code for loading the data and instructions. This is a
high-dimensional dataset with the number of samples
![n=71](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%3D71 "n=71")
much less than the number of predictors
![p=4088](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%3D4088 "p=4088").
The goal is to predict the riboflavin production rate based on gene
expression.

You will be asked to do the following:

-   use **fitLASSO** function to see how the sparsity changes with
    ![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
    value, and test the speed
-   use **cvLASSO** function to select the tuning parameter, see how
    ![CV(\lambda)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;CV%28%5Clambda%29 "CV(\lambda)")
    changes with lambda
    
    
    
## Part 2 - Reimplement in C++

In this part, I practiced Rcpp and the use of Armadillo library by implementing the LASSO coordinate-descent algorithm from Part 1. I only focus on `fitLASSOstandardized` functions and corresponding helpers in an attempt to speed up the function. 


## Basic Overview of the functions

**LassoInC.cpp** contains the starter for C++ code. You will have to modify this starter to create the following functions:

  - `soft_c` - this is the C++ version of the `soft` function from Part 1. It should return the soft-thresholded value of a. 
  
  - `lasso_c` - this is the C++ version of the `lasso` function from Part 1. It should return the value of lasso objective function at given arguments.
  
  - `fitLASSOstandardized_c` - this is a *slightly modified* C++ version of the `fitLASSOstandardized` function from Part 1. The  modification: it only returns $\beta$ vector (it no longer returns objective function value)
  
  - `fitLASSOstandardized_seq_c` - this is a *slightly modified* C++ version of the `fitLASSOstandardized` function from Part 1. The modification: it only takes used-supplied lambda sequence and it only returns matrix of $\beta$ values. You can assume that the user-supplied sequence has already been sorted from largest to smallest.
  
