---
title: spatial confounding in Cox process models
author: Wesley Brooks
---

# Spatial confounding
Spatial confounding in regression models leads to biased estimation of the model parameters. Spatial confounding may occur when the spatial structure of observations is correlated with covariates in a regression model.

# Cox process
A Cox process is a kind of point process that includes a spatial random effect. Perhaps the simplest possible example: take a homogeneous Poisson process, and then add to the log-intensity a Gaussian random field.

# Estimating regression coefficients in a point process model
Here there are many options - see the paper by Renner et al. (2015). For estimating a Cox process, not all methods will work. Think of the Cox process as a hierarchical model where the response is above the random effect in the hierarchy. This kind of model is often estimated in a Bayesian way, which can be difficult.

One goal is to estimate the Cox process with more standard software, such as lme4. Casting the Cox process as a GLM is not difficult: generate many (e.g. hundreds of thousands) of "quadrature points" that are not at point locations. Then the DWPR approach can use the quadrature points as background (i.e. assume a zero count was observed there) along with the point locations with large counts. The quadrature points are given large weight. Then use typical GLM machinery with a Poisson response to estimate the regression coefficients. 

## Estimating the regression coefficients in a Cox model
For a Cox model, can the GLM machinery be used?
- Need to cast the spatial effect as a random effect, presumable
- Use a geostatistical covariance model?

# Notes
Quadrature seems to be tiling the region and integrating the intensity on each tile. Need to multiply intensity by the area of the tile.

