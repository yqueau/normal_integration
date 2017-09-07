# Codes_Integration
Matlab codes for integration of normals over a non-rectangular 2D grid, without boundary condition.

## Introduction

In many computer vision applications (e.g. photometric stereo, shape-from-shading, shape-from-polarization or deflectometry), one estimates the local surface orientation (i.e., normals) in each pixel. A subsequent step consists in integrating these into a depth map. These Matlab codes implement the variational normal integration methods discussed in [1]. This journal paper summarizes research previously presented in the conference papers [2,3,4].  

Features:
- A fast (almost  n log(n)) quadratic integration over a non-rectangular domain, without boundary condition, without parameter to tune.
- Various non-quadratic, discontinuity-preserving integrators (slower than quadratic and a parameter needs to be tuned).
- Possibility to include a depth prior in each method (see Sec. 2.4 in [1])  

## References

[1] "Variational Methods for Normal Integration", Quéau et al., Submitted to JMIV. 2017.
[2] "Integration of a Normal Field without Boundary Condition", Durou and Courteille, ICCVW 2007
[3] "Integrating the Normal Field of a Surface in the Presence of Discontinuities", Durou et al., EMMCVPR2009
[4] "Edge-Preserving Integration of a Normal Field: Weighted Least Squares and L1 Approaches", Quéau and Durou, SSVM2015 

Please cite [1] if using the provided codes for your own research. The quadratic method is an extension of the method in [2]. The non-convex integrator was introduced in [3], and the TV one in [4]. 

Author of codes: Yvain Quéau, Technical University Munich, yvain.queau@tum.de


## Contents

The main fuctions are in the Toolbox folder:
- `make_gradient.m`: given a 2D binary mask, returns the matrix differentiation operators in all 4 directions (D_{u/v}^{+/-} in Sect. 3 in [1])  
- `smooth_integration.m`: function for quadratic integration over a non-rectangular grid (Sec. 3 in [1])
- `tv_integration.m`: function for TV integration over a non-rectangular grid (Sec. 4.1 in [1])
- `phi1_integration.m`: function for non-convex (Phi_1 estimator) integration over a non-rectangular grid (Sec. 4.2 in [1])
- `phi2_integration.m`: function for non-convex (Phi_2 estimator) integration over a non-rectangular grid (Sec. 4.2 in [1])
- `anisotropic_diffusion_integration.m`: function for anisotropic diffusion integration over a non-rectangular grid (Sec. 4.3 in [1])
- `mumford_shah_integration.m`: function for Mulford-Sjaj integration over a non-rectangular grid (Sec. 4.4 in [1])

## Usage
- All methods require to provide:
 * p: estimation of the gradient in bottom direction (matrix)
 * q: estimation of the gradient in right direction (matrix)
- Optional parameters common to all methods
 * mask: binary mask of the area of interest (matrix)
 * lambda: field of regularization weights for depth prior (matrix)
 * z0: depth prior (matrix)
- Discontinuity-preserving methods require a few other settings, see demo 4 for details. 

## Dependencies

We strongly recommend to use the CMG preconditioner from Koutis et al., which can be downloaded here: 
http://www.cs.cmu.edu/~jkoutis/cmg.html

If CMG it is not installed, set the "precond" parameter to "none". This will be (a lot) slower, but it should run without any additional library.

## Demos

The following demo files are provided: 

- `demo_1_domain.m` : demo of fast quadratic integration over a non-rectangular grid. Script shows that when explicitly using the domain, integration is faster and way more accurate (no bias on the boundary of the object due to discontinuity).

- `demo_2_CPU_eval.m` : code for comparison of asymptotic complexity of the quadratic approach with DCT and Sylvester approaches. 

- `demo_3_noise_eval.m` : code for comparison of accuracy of the quadratic approach with DCT and Sylvester approaches. 

- `demo_4_discontinuities.m` : demo of the four discontinuity-preserving methods. 


