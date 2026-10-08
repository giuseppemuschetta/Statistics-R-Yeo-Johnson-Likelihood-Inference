**Author:** Giuseppe Muschetta  
**University:** University of Pisa  
**Degree:** M.Sc. in Data Science & Business Informatics  
**Academic year:** 2024–2025

# Yeo–Johnson: From Theory to Statistical Inference

**Statistical computing for Data Science: custom power transformations, numerical optimization and likelihood-based analysis in R.**

This project connects mathematical transformation theory with practical statistical inference. Inspired by Yeo and Johnson's *A New Family of Power Transformations to Improve Normality or Symmetry*, it implements power transformations directly in R, investigates their behaviour on synthetic distributions, and develops an inferential analysis of Darwin's classic plant-growth experiment.

The work brings together **distribution geometry, parameter estimation, assumption diagnostics and hypothesis testing** in a coherent investigation of how transformations affect statistical modelling.

## Project at a glance

| Dimension | Scope |
| :--- | :--- |
| Statistical foundation | Yeo–Johnson, Box–Cox and signed-power transformations |
| Synthetic experiments | 10,000 sampled observations per mixture example, drawn from a 3,000-point grid |
| Numerical methods | Parameter selection with L-BFGS-B and joint likelihood optimization |
| Applied case study | 15 paired plant-height differences from Darwin's fertilization experiment |
| Inferential methods | One-sided t-tests, Shapiro–Wilk diagnostics, skewness, likelihood-ratio analysis and Hessian calculations |
| Implementation | Four R scripts, with shared statistical functions and comparative visualizations |

## The statistical motivation

Transformations can help reduce skewness and make distributional assumptions more plausible. Their usefulness depends on the shape of the data, the domain of the transformation and the statistical question being asked.

Box–Cox is defined for strictly positive observations. Yeo–Johnson accommodates **negative values, zero and positive values** within a single parameterized family. This makes it particularly relevant to Data Science problems involving signed measurements, differences between paired observations or features whose distributions depart from symmetry.

The project studies these choices through three connected stages.

## 1. Investigating the limits of power transformations

[`1st_part.R`](1st_part.R) introduces a deliberately asymmetric distribution: a mixture of **30% standard normal density** and **70% Gamma(shape = 4, rate = 1), shifted by −2**.

The script samples 10,000 observations from a discretized density, selects a signed-power transformation parameter through bounded numerical optimization, and compares:

- The original mixture density.
- A kernel-density estimate of the transformed observations.
- A Gaussian density fitted using the transformed sample's mean and standard deviation.

This stage examines the relationship between a transformation's shape and the distribution it produces, establishing the motivation for the Yeo–Johnson family.

## 2. Exploring Yeo–Johnson geometry and distribution shaping

[`2nd_part.R`](2nd_part.R) studies the transformation across **λ = 0, 0.5, 1, 1.5 and 2**, compares it with Box–Cox, and constructs visualizations inspired by the original paper.

The analysis covers the identity case at λ = 1, the different curvature for λ below and above 1, and the behaviour of the positive and negative branches. A second synthetic-mixture example uses numerical parameter selection to compare the transformed density with a fitted Gaussian.

The central transformation is implemented directly in [`functions.R`](functions.R), including its logarithmic special cases.

### The transformation behind the implementation

For an observation $x$, the Yeo–Johnson transformation is

```math
T_\lambda(x)=
\begin{cases}
\dfrac{(x+1)^\lambda-1}{\lambda}, & x\geq 0,\ \lambda\neq 0,\\[6pt]
\log(x+1), & x\geq 0,\ \lambda=0,\\[6pt]
-\dfrac{(1-x)^{2-\lambda}-1}{2-\lambda}, & x<0,\ \lambda\neq 2,\\[6pt]
-\log(1-x), & x<0,\ \lambda=2.
\end{cases}
```

At $\lambda=1$, $T_\lambda(x)=x$. The same family therefore includes the untransformed observations as a reference case.

## 3. From transformed data to statistical inference

[`3rd_part.R`](3rd_part.R) applies the framework to **Darwin's experiment on cross-fertilization versus self-fertilization**, using the height differences from 15 matched plant pairs.

The analysis follows the statistical question through a complete sequence:

1. **Assess the original observations:** a one-sided t-test on the paired differences, a skewness calculation, a Shapiro–Wilk test and a density plot.
2. **Estimate transformation and distribution parameters:** numerical optimization of λ, followed by joint refinement of λ, μ and σ² using the transformed-data likelihood.
3. **Examine parameter uncertainty:** calculation of a covariance matrix from the numerical Hessian.
4. **Reassess the transformed observations:** repeat the normality diagnostics and the one-sided t-test.
5. **Compare transformation choices:** a likelihood-ratio statistic for the optimized model against λ = 1, with a χ² reference distribution.
6. **Return to the original measurement scale:** a back-transformed lower 0.01 quantile estimate based on a Student-t reference.

The likelihood includes the **Jacobian contribution**, connecting the transformation to the probability model for the original observations.

### Likelihood formulation

Under a Gaussian model for the transformed observations, the log-likelihood is

```math
\begin{aligned}
\ell(\lambda,\mu,\sigma^2\mid x)
&=-\frac{n}{2}\log(2\pi)-\frac{n}{2}\log(\sigma^2)\\[6pt]
&\quad-\frac{1}{2\sigma^2}\sum_{i=1}^{n}\left[T_\lambda(x_i)-\mu\right]^2\\[6pt]
&\quad+(\lambda-1)\sum_{i=1}^{n}\mathrm{sgn}(x_i)\log(1+|x_i|).
\end{aligned}
```

The shared `log_likelihood()` function returns the negative log-likelihood for minimization with R's optimization routines.

## Technical skills demonstrated

| Skill | Evidence in the project |
| :--- | :--- |
| Translating statistical theory into code | Direct implementation of transformation formulas and a likelihood with its Jacobian term |
| Numerical computing | Bounded parameter optimization, joint estimation and numerical Hessian calculations |
| Statistical reasoning | Connecting distributional assumptions, diagnostics and hypothesis tests |
| Scientific visualization | Comparisons of transformation curves, mixture densities and Gaussian reference densities |
| Modular R programming | Shared statistical utilities, documented functions and three focused analysis scripts |

## Explore the repository

| File | Role |
| :--- | :--- |
| [`functions.R`](functions.R) | Shared transformations, inverse and derivative helpers, numerical property diagnostics, likelihood, skewness and quantile calculations |
| [`1st_part.R`](1st_part.R) | Signed-power exploration on a synthetic asymmetric mixture |
| [`2nd_part.R`](2nd_part.R) | Yeo–Johnson geometry, Box–Cox comparison and transformed-mixture analysis |
| [`3rd_part.R`](3rd_part.R) | Darwin's paired-data case study and likelihood-based inference |
| [`paper_TransformNormality.pdf`](paper_TransformNormality.pdf) | Original research paper by Yeo and Johnson |

**Tools:** R, its built-in statistical and graphics functions, `ggplot2` and `gridExtra`. The examples use generated observations or the 15 historical height differences embedded in the source. The scripts produce plots and printed statistics within the R session.

## Reference

Yeo, I.-K., & Johnson, R. A. (2000). **A New Family of Power Transformations to Improve Normality or Symmetry.** *Biometrika*, 87(4), 954–959. [Read the original publication](https://www.jstor.org/stable/2673623).

This repository presents Giuseppe Muschetta's R implementation and statistical exploration of the published framework.
