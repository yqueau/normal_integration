function [Z,Zx,Zy,x,y] = g2sTestSurf(m,n,nodeType,plotSurf) 
%
% Purpose : Generates a test surface along with its analytic gradient field
%   for testing surface reconstruction from gradient field algorithms.
%
% Use (syntax):
%   [Z,Zx,Zy,x,y] = g2s( m, n )
%   [Z,Zx,Zy,x,y] = g2s( m, n, nodeType )
%   [Z,Zx,Zy,x,y] = g2s( m, n, nodeType, plotSurf )
%
% Input Parameters :
%   m, n := the matrix size of the sampled surface
%   nodeType := 'even'  (Evenly spaced points)
%               'cheb'  (Chebyshev-like points, scaled to include end-points)
%               'rand'  (Randomized points)
%   plotSurf := 0/1, plot the surface and its gradient or not.
%
% Return Parameters :
%   Z := The sampled test surface
%   Zx, Zy := Components of the discrete gradient field
%   x, y := support vectors of nodes of the domain of the gradient
%
% Description and algorithms:
%   The algorithm samples a test function of the form:
%   z(x,y) = sum_{k=1}^{n} alpha_k exp( -[x-X_k y-Y_k] * A [x-X_k ; y-Y_k] ) 
%   along with its analytic partial derivatives.
%
% References :
%    @inproceedings{
%    Harker2008c,
%       Author = {Harker, M. and O'Leary, P.},
%       Title = {Least Squares Surface Reconstruction from Measured Gradient Fields},
%       BookTitle = {CVPR 2008},
%       Address= {Anchorage, AK},
%       Publisher = {IEEE},
%       Pages = {1-7},
%          Year = {2008} }
%
%    @inproceedings{
%    harker2011,
%       Author = {Harker, M. and O'Leary, P.},
%       Title = {Least Squares Surface Reconstruction from Gradients:
%           \uppercase{D}irect Algebraic Methods with Spectral, \uppercase{T}ikhonov, and Constrained Regularization},
%       BookTitle = {IEEE CVPR},
%       Address= {Colorado Springs, CO},
%       Publisher = {IEEE},
%       Pages = {2529--2536},
%          Year = {2011} }
%
% Author :  Matthew Harker and Paul O'Leary
% Date :    25. June 2013
% Version : 1.0
%
% (c) 2013 Matthew Harker and Paul O'Leary, 
% Chair of Automation, University of Leoben, Leoben, Austria
% email: office@harkeroleary.org, 
% url: www.harkeroleary.org
%
% History:
%   Date:           Comment:
%   Apr. 17, 2013   Original Version
%
%================================
% Argument Check:
%================================
%
if nargin == 2
    %
    plotSurf = 0 ;
    nodeType = 'even' ;
    %
elseif nargin == 3
    %
    plotSurf = 0 ;
    %
end
%
nodeType = lower( nodeType ) ;
%
%===============================
% Set up Nodes:
%===============================
%
a = -1 ;
b = 10 ;
c = -1 ;
d = 10 ;
%
% Scale a given node-spacing onto the intervals:
%
switch nodeType
    case 'even'
        %
        x = linspace(a,b,n)' ;
        y = linspace(c,d,m)' ;
        %
    case 'cheb'
        %
        nodesX = - cos( pi * linspace(0,1,m) );
        nodesY = - cos( pi * linspace(0,1,n) );
        x = a + (b-a)*( nodesX + 1 )/2 ;
        y = c + (d-a)*( nodesY + 1 )/2 ;
        %
    case 'rand'
        %
        x = linspace(a,b,n)' ;
        y = linspace(a,b,m)' ;
        hx = (b-a)/(n-1) ;
        hy = (b-a)/(m-1) ;
        x(2:n-1) = x(2:n-1) + (hx/4)*rand(n-2,1) ;
        y(2:m-1) = y(2:m-1) + (hy/4)*rand(m-2,1) ;
        %
    otherwise
        %
        error('Not a valid node type') ;
        %
end
%
%=================================================
% Generate the Analytic Surface and its Gradient
%=================================================
%
[X,Y] = meshgrid( x, y ) ;
%
Z = 0.5e1 / 0.2e1 * exp(-0.3e1 / 0.16e2 * X .^ 2 + 0.5e1 / 0.8e1 * X - 0.19e2 / 0.16e2 - Y .* X / 0.8e1 + 0.7e1 / 0.8e1 * Y - 0.3e1 / 0.16e2 * Y .^ 2) + 0.3e1 * exp(-0.2e1 / 0.7e1 * X .^ 2 + 0.32e2 / 0.7e1 * X - 0.142e3 / 0.7e1 - Y .* X / 0.7e1 + 0.15e2 / 0.7e1 * Y - Y .^ 2 / 0.7e1) - 0.5e1 * exp(-0.5e1 / 0.18e2 * X .^ 2 + 0.20e2 / 0.9e1 * X - 0.125e3 / 0.18e2 + Y .* X / 0.9e1 + 0.5e1 / 0.9e1 * Y - Y .^ 2 / 0.9e1) - 0.2e1 * exp(-0.3e1 / 0.28e2 * X .^ 2 - X / 0.7e1 - 0.75e2 / 0.7e1 + Y .* X / 0.14e2 + 0.19e2 / 0.7e1 * Y - 0.5e1 / 0.28e2 * Y .^ 2) + 0.5e1 * exp(-X .^ 2 / 0.6e1 + 0.14e2 / 0.3e1 * X - 0.194e3 / 0.3e1 - Y .* X / 0.3e1 + 0.38e2 / 0.3e1 * Y - 0.2e1 / 0.3e1 * Y .^ 2) ;
Zx = 0.5e1 / 0.2e1 * (-0.3e1 / 0.8e1 * X + 0.5e1 / 0.8e1 - Y / 0.8e1) .* exp(-0.3e1 / 0.16e2 * X .^ 2 + 0.5e1 / 0.8e1 * X - 0.19e2 / 0.16e2 - Y .* X / 0.8e1 + 0.7e1 / 0.8e1 * Y - 0.3e1 / 0.16e2 * Y .^ 2) + 0.3e1 * (-0.4e1 / 0.7e1 * X + 0.32e2 / 0.7e1 - Y / 0.7e1) .* exp(-0.2e1 / 0.7e1 * X .^ 2 + 0.32e2 / 0.7e1 * X - 0.142e3 / 0.7e1 - Y .* X / 0.7e1 + 0.15e2 / 0.7e1 * Y - Y .^ 2 / 0.7e1) - 0.5e1 * (-0.5e1 / 0.9e1 * X + 0.20e2 / 0.9e1 + Y / 0.9e1) .* exp(-0.5e1 / 0.18e2 * X .^ 2 + 0.20e2 / 0.9e1 * X - 0.125e3 / 0.18e2 + Y .* X / 0.9e1 + 0.5e1 / 0.9e1 * Y - Y .^ 2 / 0.9e1) - 0.2e1 * (-0.3e1 / 0.14e2 * X - 0.1e1 / 0.7e1 + Y / 0.14e2) .* exp(-0.3e1 / 0.28e2 * X .^ 2 - X / 0.7e1 - 0.75e2 / 0.7e1 + Y .* X / 0.14e2 + 0.19e2 / 0.7e1 * Y - 0.5e1 / 0.28e2 * Y .^ 2) + 0.5e1 * (-X / 0.3e1 + 0.14e2 / 0.3e1 - Y / 0.3e1) .* exp(-X .^ 2 / 0.6e1 + 0.14e2 / 0.3e1 * X - 0.194e3 / 0.3e1 - Y .* X / 0.3e1 + 0.38e2 / 0.3e1 * Y - 0.2e1 / 0.3e1 * Y .^ 2);
Zy = 0.5e1 / 0.2e1 * (-X / 0.8e1 + 0.7e1 / 0.8e1 - 0.3e1 / 0.8e1 * Y) .* exp(-0.3e1 / 0.16e2 * X .^ 2 + 0.5e1 / 0.8e1 * X - 0.19e2 / 0.16e2 - Y .* X / 0.8e1 + 0.7e1 / 0.8e1 * Y - 0.3e1 / 0.16e2 * Y .^ 2) + 0.3e1 * (-X / 0.7e1 + 0.15e2 / 0.7e1 - 0.2e1 / 0.7e1 * Y) .* exp(-0.2e1 / 0.7e1 * X .^ 2 + 0.32e2 / 0.7e1 * X - 0.142e3 / 0.7e1 - Y .* X / 0.7e1 + 0.15e2 / 0.7e1 * Y - Y .^ 2 / 0.7e1) - 0.5e1 * (X / 0.9e1 + 0.5e1 / 0.9e1 - 0.2e1 / 0.9e1 * Y) .* exp(-0.5e1 / 0.18e2 * X .^ 2 + 0.20e2 / 0.9e1 * X - 0.125e3 / 0.18e2 + Y .* X / 0.9e1 + 0.5e1 / 0.9e1 * Y - Y .^ 2 / 0.9e1) - 0.2e1 * (X / 0.14e2 + 0.19e2 / 0.7e1 - 0.5e1 / 0.14e2 * Y) .* exp(-0.3e1 / 0.28e2 * X .^ 2 - X / 0.7e1 - 0.75e2 / 0.7e1 + Y .* X / 0.14e2 + 0.19e2 / 0.7e1 * Y - 0.5e1 / 0.28e2 * Y .^ 2) + 0.5e1 * (-X / 0.3e1 + 0.38e2 / 0.3e1 - 0.4e1 / 0.3e1 * Y) .* exp(-X .^ 2 / 0.6e1 + 0.14e2 / 0.3e1 * X - 0.194e3 / 0.3e1 - Y .* X / 0.3e1 + 0.38e2 / 0.3e1 * Y - 0.2e1 / 0.3e1 * Y .^ 2);
%
Z = Z - mean( Z(:) ) ;
%
if plotSurf == 1
    %
    % Plot the surface and its gradient field:
    %
    subplot(1,3,1)
    surfl(x,y,Z);
    shading interp
    colormap(gray);
    axis equal
    title('Surface') ;
    xlabel('{\it x}') ;
    ylabel('{\it y}') ;
    subplot(1,3,2)
    surfl(x,y,Zx);
    shading interp
    colormap(gray);
    axis equal
    title('{\it x}-Derivative') ;
    xlabel('{\it x}') ;
    ylabel('{\it y}') ;
    subplot(1,3,3)
    surfl(x,y,Zy);
    shading interp
    colormap(gray);
    axis equal
    title('{\it y}-Derivative') ;
    xlabel('{\it x}') ;
    ylabel('{\it y}') ;
    %
end
%
%========
% END
%========
