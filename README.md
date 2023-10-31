
# Homework 4+5: Coordinate-descent algorithm for LASSO

## Introduction

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
respectively. We would like to fit linear model

![Y = \beta_0 + X\beta + \varepsilon,](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%20%3D%20%5Cbeta_0%20%2B%20X%5Cbeta%20%2B%20%5Cvarepsilon%2C "Y = \beta_0 + X\beta + \varepsilon,")

where the sample size
![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
is small compared to the number of covariates
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p").
We will use LASSO algorithm to fit this model (find
![\beta_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_0 "\beta_0")
and
![\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")),
and use 5-fold cross-validation to select the tuning parameter.


## Part 1 - Implement LASSO with cross validation

For this part you need to implement LASSO and 5 fold cross validation to select the tuning parameter in R.

1)  We will center
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
    (because of scaling). See class notes for more.

2)  We will solve the following LASSO problem for various
    ![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
    values

    ![\widetilde \beta = \arg\min\_{\beta}\left\\{(2n)^{-1}\\\|\widetilde Y-\widetilde X\beta\\\|\_2^2 + \lambda \\\|\beta\\\|\_1\right\\}.](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cwidetilde%20%5Cbeta%20%3D%20%5Carg%5Cmin_%7B%5Cbeta%7D%5Cleft%5C%7B%282n%29%5E%7B-1%7D%5C%7C%5Cwidetilde%20Y-%5Cwidetilde%20X%5Cbeta%5C%7C_2%5E2%20%2B%20%5Clambda%20%5C%7C%5Cbeta%5C%7C_1%5Cright%5C%7D. "\widetilde \beta = \arg\min_{\beta}\left\{(2n)^{-1}\|\widetilde Y-\widetilde X\beta\|_2^2 + \lambda \|\beta\|_1\right\}.")

    To solve LASSO, we will use coordinate-descent algorithm with **warm
    starts**.

3)  We will use the
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
    based on back-scaling and back-centering (see class notes).
    

## Starter code

The starter code for all functions with detailed description is provided
in **LassoFunctions.R**. The described functions should have
input/output **exactly as specified**, but you are welcome to create any
additional functions you need. You are not allowed to use any outside
libraries for these functions. I encourage you to work gradually through
the functions and perform frequent testing as many functions rely on the
previous ones.

Things to keep in mind when implementing:

-   You want to make sure that parameters supplied to one function are
    correctly used in subsequent functions (i.e. the convergence level
    ![\varepsilon](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cvarepsilon "\varepsilon"))

-   You should check your code on simple examples before proceeding to
    the data (i.e. what happens when large lambda is supplied? What
    happens on toy example used in class?). I will use automatic tests
    to check that your code is correct on more than just the data
    example with different combinations of parameters.

-   I will test your functions speed on large
    ![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p"),
    small
    ![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
    dataset (below), so you need to use coordinate-descent
    implementation that is optimized for this scenario to pass the speed
    requirements (see class notes). I recommend writing a correct code
    first using whichever version of algorithm you find the easiest, and
    then optimizing for speed (this will help avoid mistakes).

-   I expect it will take you some time to figure out how to split the
    data for cross-validation. Keep in mind that the split should be
    random, in roughly equal parts, and should work correctly with any
    sample size
    ![n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n "n")
    and any integer number of folds
    ![K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K "K")
    as long as
    ![n\geq K](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%5Cgeq%20K "n\geq K").

-   `glmnet` function within glmnet R package is probably the most
    popular LASSO solver, however the outputs are not directly
    comparable as glmnet does additional scaling of
    ![Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y "Y")
    to have unit variance with 1/n formula, that is
    ![\tilde Y](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctilde%20Y "\tilde Y")
    for glmnet will satisfy
    ![n^{-1}\tilde Y^{\top}\tilde Y = 1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%5E%7B-1%7D%5Ctilde%20Y%5E%7B%5Ctop%7D%5Ctilde%20Y%20%3D%201 "n^{-1}\tilde Y^{\top}\tilde Y = 1")
    (see Details in ? `glmnet`). However, if you impose the same
    standardization, you may be able to use glmnet for testing. In terms
    of speed, **glmnet** is very highly optimized with underlying
    Fortran so it would be very fast compared to your code. On my
    computer, glmnet with 60 tuning parameters on riboflavin dataset
    takes around 42 milliseconds.

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

In this part, you will practice Rcpp and the use of Armadillo library by implementing the LASSO coordinate-descent algorithm from Part 1. We will only focus on `fitLASSOstandardized` functions and corresponding helpers as you will see the biggest speed improvement there. For simplicity, we will also **avoid doing any compatibility/input checks** in C++. In practice, you will often have an R code wrapper with all compatibility checks that then calls the corresponding function in C++.


## Starter code

**LassoInC.cpp** contains the starter for C++ code. You will have to modify this starter to create the following functions:

  - `soft_c` - this is the C++ version of the `soft` function from Part 1. It should return the soft-thresholded value of a. 
  
  - `lasso_c` - this is the C++ version of the `lasso` function from Part 1. It should return the value of lasso objective function at given arguments.
  
  - `fitLASSOstandardized_c` - this is a *slightly modified* C++ version of the `fitLASSOstandardized` function from Part 1. The  modification: it only returns $\beta$ vector (it no longer returns objective function value)
  
  - `fitLASSOstandardized_seq_c` - this is a *slightly modified* C++ version of the `fitLASSOstandardized` function from Part 1. The modification: it only takes used-supplied lambda sequence and it only returns matrix of $\beta$ values. You can assume that the user-supplied sequence has already been sorted from largest to smallest.
  
All these functions should be directly available from R after executing **sourceCpp** command with corresponding source file (see **TestsCpp.R**). You are welcome to create additional C++ functions within the same file if you want, we will only explicitly test the above 4 functions. You can assume that all the inputs to these functions are correct (e.g. dimensions match, no negative tuning parameters, etc.)

**TestsCpp.R**  contains starter file for testing your C++ functions against your R functions. You should do the following

  - develop at least 2 tests for each function **to check equality of the corresponding outputs between your R and C++ versions**
  - do 2 microbenchmarks comparisons as per comments
  - do speed test on riboflavin data as coded in the end of that file
  
Things to keep in mind when implementing:

 - you will use Armadillo library classes for matrices and vectors, and because of this, the whole Armadillo library of functions is at your disposal. As it's impossible to cover all functions in class, a part of this assignment is learning how to search for needed functions on your own using [Armadillo documentation](http://arma.sourceforge.net/docs.html). Using Ctrl(Command) + F on that page with a key word for the function you need, e.g. "inverse", will help you narrow down the function list for what you need.
 
 - one of the best ways to quickly learn how to translate the code is to study some relevant examples. In addition to examples in class, you may find the following [gallery of Rcpp and RcppArmadillo examples](https://gallery.rcpp.org) useful

 - test small chunks before testing larger chunks as it's harder to debug C++ code, and it's easier to miss a mistake in C++. I recommend first testing that a particular line of code does exactly what you think it does before going to the whole function
 
 - there are many more ways in C++ code to achieve a similar outcome than in R so we will be expecting to see much more differences across the assignments. If your code is heavily based on particular examples you found online, you will be expected to indicate this in the comments


## Grading for this assignment

Your assignment will be graded based on

-   correctness *(50% of the grade)*

Take advantage of objective function values over iterations as a way to
indirectly check the correctness of your function. Also recall that you
know the right solution in special cases, so you can check your function
in those cases (i.e. when
![\lambda](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda "\lambda")
is very large, or when
![\lambda = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda%20%3D%200 "\lambda = 0")
and you have a nice problem)

-   speed of implementations (use warm starts and coordinate descent
    optimized for large p settings) *(30% of the grade)*

You will get full points if your code is **at most twice slower**
comparable to mine (fitLASSO with 60 tuning parameters on riboflavin
data takes around 5.8 seconds on my laptop). You will loose 5 points for
every fold over. You will get +5 bonus points if your **completely
correct** code on 1st submission is faster than mine (median your
time/median mine time \< 0.9).

-   code style/documentation *(10% of the grade)*

You need to comment different parts of the code so it’s clear what they
do, have good indentation, readable code with names that make sense. See
guidelines on R style, and posted grading rubric.

-   version control/commit practices *(10% of the grade)*

I expect you to start early on this assignment, and work gradually. You
want to commit often, have logically organized commits with short
description that makes sense. See guidelines on good commit practices,
and posted grading rubric.
