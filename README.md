# Codes_Integration
Matlab codes for integration of normals (gradient) over a non-rectangular 2D grid, without boundary condition.

## Introduction

In many computer vision applications (e.g. photometric stereo, shape-from-shading, shape-from-polarization or deflectometry), one estimates the local surface orientation (i.e., normals) in each pixel. A subsequent step consists in integrating these normals into a depth map. These Matlab codes implement the variational normal integration methods discussed in [1], and three famous methods for solving the Poisson equation, which were discussed in our survey [2]. These journal papers summarize research previously presented in the conference papers [3,4,5].  

Features:
- Implementations of Poisson solvers using FFT (see Sec. 3.3 in [2]), using DCT (see Sec. 3.4 in [2], and the modified Jacobi scheme of Horn and Brooks (see Sec. 3.2 in [2], and demo 1);   
- A fast (almost  n log(n)) quadratic integration over a non-rectangular domain, without boundary condition, without parameter to tune (see Sec. 3 in [1], and demo 2);
- Various non-quadratic, discontinuity-preserving integrators: total variation, non-convex, anisotropic diffusion and Mumford-Shah (all are slower than quadratic integration, and require at least one parameter to be tuned) (see Sec. 4 in [1], and demo 3);
- Possibility to include a depth prior in each method (see Sec. 2.4 in [1]).


## Demos

The following demo files are provided: 

- `demo_1_survey.m` : demo of the methods presented in Sec. 3 of the survey paper [2]. This shows the importance of boundary conditions, hence the superiority of DCT over DST/FFT for solving the Poisson equation. However, it handles only a rectangular domain, hence Horn and Brook's method is much more accurate for non-rectangular domains. Still, the latter is very slow, hence the quadratic method proposed in [1] is much better: it handles free boundary and free-form domain, while being almost as fast as DCT.   

- `demo_2_domain.m` : demo of fast quadratic integration over a non-rectangular grid. Script shows that when explicitly using the domain, integration is faster and way more accurate (no bias on the boundary of the object due to discontinuity).

- `demo_3_discontinuities.m` : demo of the four discontinuity-preserving methods. They can be used if the domain of integration has not been pre-calculated. 



## Contents

The main fuctions for the new variational methods in [1] are in the Toolbox/ folder:
- `make_gradient.m`: given a 2D binary mask, returns the matrix differentiation operators in all 4 directions (D_{u/v}^{+/-} in Sec. 3 in [1])  
- `smooth_integration.m`: function for quadratic integration over a non-rectangular grid (Sec. 3 in [1])
- `tv_integration.m`: function for TV integration over a non-rectangular grid (Sec. 4.1 in [1])
- `phi1_integration.m`: function for non-convex (Phi_1 estimator) integration over a non-rectangular grid (Sec. 4.2 in [1])
- `phi2_integration.m`: function for non-convex (Phi_2 estimator) integration over a non-rectangular grid (Sec. 4.2 in [1])
- `anisotropic_diffusion_integration.m`: function for anisotropic diffusion integration over a non-rectangular grid (Sec. 4.3 in [1])
- `mumford_shah_integration.m`: function for Mulford-Sjaj integration over a non-rectangular grid (Sec. 4.4 in [1])

The four Poisson solvers discussed in [2, Sec. 3] are also provided:
- `horn_brooks.m`: implementation of the modified Horn and Brook's scheme (Jacobi iterations) for Poisson integration over a non-rectangular grid. Needs no boundary condition, but very slow (Sec. 3.2 in [2]).
- `FFT_Poisson.m`: implementation of the FFT integrator. Super fast, but requires a rectangular grid and periodic boundary condition (Sec. 3.3 in [2])
- `DCT_Poisson`: implementation of the DCT integrator. Still very fast, and requires no boundary condition, but domain must be rectangular (Sec. 3.4 in [2])
- `DST_Poisson`: implementation of the DST integrator. Still very fast, handles Dirichlet boundary condition. Domain must be rectangular (Sec. 3.4 in [2]). 


## Dependencies

We strongly recommend to use the CMG preconditioner from Koutis et al., which can be downloaded here: 
http://www.cs.cmu.edu/~jkoutis/cmg.html

If CMG it is not installed, set the "precond" parameter to "none" (no preconditioning, can be very slow for large data) or to "ichol" (modified incomplete Cholesky preconditioner advised in [6], slower than CMG but overall OK). This will be slower, but it should run without any additional library.



## Usage
- All methods require to provide:
 * p: estimation of the gradient in bottom direction (matrix)
 * q: estimation of the gradient in right direction (matrix)
- Optional parameters common to all methods
 * mask: binary mask of the area of interest (matrix)
 * lambda: field of regularization weights for depth prior (matrix)
 * z0: depth prior (matrix)
- Discontinuity-preserving methods and Horn and Brook's one require a few other settings, see demos 1 and 3 for details. 

## References

[1] "Variational Methods for Normal Integration", Quéau et al., Journal of Mathematical Imaging and Vision 60(4), pp 609--632, 2018. (Arxiv preprint: https://arxiv.org/abs/1709.05965)

[2] "Normal Integration: a Survey", Quéau et al., Journal of Mathematical Imaging and Vision 60(4), pp 576--593, 2018. (Arxiv preprint: https://arxiv.org/abs/1709.05940)

These methods build upon three previous conference papers. The new quadratic method is an extension of the method in [3]. The non-convex integrator was introduced in [4], and the TV one in [5]. 

[3] "Integration of a Normal Field without Boundary Condition", Durou and Courteille, ICCVW 2007

[4] "Integrating the Normal Field of a Surface in the Presence of Discontinuities", Durou et al., EMMCVPR2009

[5] "Edge-Preserving Integration of a Normal Field: Weighted Least Squares and L1 Approaches", Quéau and Durou, SSVM2015 

The modified incomplete Cholesky preconditioner is advised in: 

[6] "Fast and accurate surface normal integration on non-rectangular domains", Bähr et al., Computational Visual Media 3(2), pp. 107--129, 2017


Author of codes: Yvain Quéau, Technical University Munich, yvain.queau@tum.de



