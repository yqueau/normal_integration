function Z = harker_oleary( Zy, Zx, x, y, N )
%
% Purpose : Computes the Global Least Squares reconstruction of a surface
%   from its gradient field.
%
% Use (syntax):
%   Z = g2s( Zx, Zy, x, y )
%   Z = g2s( Zx, Zy, x, y, N )
%
% Input Parameters :
%   Zx, Zy := Components of the discrete gradient field
%   x, y := support vectors of nodes of the domain of the gradient
%   N := number of points for derivative formulas (default=3)
%
% Return Parameters :
%   Z := The reconstructed surface
%
% Description and algorithms:
%   The algorithm solves the normal equations of the Least Squares cost
%   function, formulated by matrix algebra:
%   e(Z) = || D_y * Z - Zy ||_F^2 + || Z * Dx' - Zx ||_F^2
%   The normal equations are a rank deficient Sylvester equation which is
%   solved by means of Householder reflections and the Bartels-Stewart
%   algorithm.
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
% Date :    17. January 2013
% Version : 1.0
%
% (c) 2013 Matthew Harker and Paul O'Leary, 
% Chair of Automation, University of Leoben, Leoben, Austria
% email: office@harkeroleary.org, 
% url: www.harkeroleary.org
%
% History:
%   Date:           Comment:
%   Feb. 9, 2011    Original Version
%

if ~all( size(Zx)==size(Zy) )
    %
    error('Gradient components must be the same size') ;
    %
end
%

if(nargin==2)
	x = transpose(1:size(Zx,2));
	y = transpose(1:size(Zx,1));	
	N = 3;
end

if ~( size(Zx,2)==length(x) ) || ~( size(Zx,1)==length(y) )
    %
    error('Support vectors must have the same size as the gradient') ;
    %
end
%


if nargin==4
    %
    N = 3 ;
    %
end
%
[m,n] = size( Zx ) ;
%
Dx = dopDiffLocal( x, N, N, 'sparse' ) ;
Dy = dopDiffLocal( y, N, N, 'sparse' ) ;
%
Z = g2sSylvester( Dy, Dx, Zy, Zx, ones(m,1), ones(n,1) ) ;
%
%========
% END
%========
%
end

function Phi = g2sSylvester( A, B, F, G, u, v )
%
% Purpose : Solves the semi-definite Sylvester Equation of the form
%   A'*A * Phi + Phi * B'*B - A'*F - G*B = 0,
%   Where the null vectors of A and B are known to be
%   A * u = 0
%   B * v = 0
%
% Use (syntax):
%   Phi = g2sSylvester( A, B, F, G, u, v )
%
% Input Parameters :
%   A, B, F, G := Coefficient matrices of the Sylvester Equation
%   u, v := Respective null vectors of A and B
%
% Return Parameters :
%   Phi := The minimal norm solution to the Sylvester Equation
%
% Description and algorithms:
%   The rank deficient Sylvester equation is solved by means of Householder
%   reflections and the Bartels-Stewart algorithm.  It uses the MATLAB
%   function "lyap", in reference to Lyapunov Equations, a special case of
%   the Sylvester Equation.
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
%   Feb. 9, 2011    Original Version
%
%----------------------
% Householder Vectors:
%----------------------
%
m = length(u) ;
n = length(v) ;
%
u(1) = u(1) + norm( u, 2 ) ;
u = ( sqrt(2) / norm( u, 2 ) ) * u ;
%
v(1) = v(1) + norm( v, 2 ) ;
v = ( sqrt(2) / norm( v, 2 ) ) * v ;
%
%--------------------------------
% Apply the Householder Updates:
%--------------------------------
%
% With:
% Pa = eye(m) - u * u' ;
% Pb = eye(n) - v * v' ;
% 
% Compute:
% A = A * Pa ;
% B = B * Pb ;
% F = F * Pb ;
% G = Pa' * G ;
%
A = A - ( A * u ) * u' ;
B = B - ( B * v ) * v' ;
F = F - ( F * v ) * v' ;
G = G - u * ( u' * G ) ;
%
%--------------------------------
% Solve the System of Equations:
%--------------------------------
%
Phi = zeros(m,n) ;
%
Phi(1,2:n) = G(1,:) / B(:,2:n)' ;
Phi(2:m,1) = A(:,2:m) \ F(:,1) ;
Phi(2:m,2:n) = lyap( A(:,2:m)'*A(:,2:m), B(:,2:n)'*B(:,2:n), -A(:,2:m)'*F(:,2:end) - G(2:end,:)*B(:,2:n) ) ;
%
%---------------------------------
% Invert the Householder Updates:
%---------------------------------
%
Phi = Phi - u * ( u' * Phi ) ;
Phi = Phi - ( Phi * v ) * v' ;
%
%Phi = Phi - u * ( u' * Phi ) - ( Phi * v ) * v' + u * ( u' * Phi * v ) * v' ;
%
%=====
% END
%=====
end

function S = dopDiffLocal( x, ls, noBfs, option )
%
% Purpose : This function generates a global matrix operator which
% implements the computation of local differentials where the vector of x
% values may be irregullarly spaced.
%
% In general the support length ls should be an odd number. There is an
% exceltion made upto to ls = 20 and ls = noPoints, in this case a full
% differentiating matrix is computed.
%
% Use (syntax):
%   S = dopDiffLocal( x, ls, noBfs, option )
%
% Input Parameters :
%   x : The vestor of x value fo rthe computation.
%   ls : The support length used for the local differential
%   noBfs : the number of basis functions to be used.
%   option: 'sparse' generated sparse matrix notations, default is full.
%
% Return Parameters :
%   S: The local differential operator
%
% Description and algorithms:
%   Local discrete orthogonal polynomials are used to generate the local approximations
%   for the dreivatives
%
%
% Author :  Matthew Harker and Paul O'Leary
% Date :    17. January 2012
% Version : 1.0
%
% (c) 2013 Matthew Harker and Paul O'Leary,
% Chair of Automation, University of Leoben, Leoben, Austria
% email: office@harkeroleary.org,
% url: www.harkeroleary.org
%
% History:
%   Date:           Comment:
%

%-----------------------------------------------------------------------
[noPts, mt] = size(x);
% Test the input paramaters
%----------------------------------------------------------------
% Use sparse matrices if necessary
%
if nargin == 4
    intOption = option;
    genSparse = true;
else
    intOption = 'full';
    genSparse = false;
end;
%
%
if mt > 1
    error('A column vector is expected for x');
end;
%
%
% Test the degree and support length for campatability
%
if noBfs > ls
    error('The number of basis functions must fulfill n <= ls');
end;
%
if ls > 13
    warning('With a support length greater than 13 there may be problems with the Runge phenomena.');
end;
%
% Compute a full matrix
%
if (ls == noPts)
    [Gt, dGt] = dop( x, noBfs );
    S = dGt * Gt';
    rS = rank( S );
    if rS < noPts - 1
        warning(['The rank of S is ',int2str(rS),' while x has n = ',int2str(noPts),' points.']);
    end;
    return;
end;
%
% Test if the support length is compatible with the number of points
% requested.
%
if noPts < ls
    error('The number of nodes n must be greater that the support length ls');
end;
%
if isEven( ls )
    error('this function is only implemented for even values of ls.');
else
    %
    %------------------------------------------------------------------------
    %
    rows = [];
    cols = [];
    vals = [];
    %
    % Determine the half length of ls this determine the upper ane lower
    % postions of Si.
    %
    ls2 = round( (ls + 1 )/2 );
    %
    % generatethe top of Si
    %
    range = (1:ls)';
    halfRange = (1:ls2)';
    
    startX = x(range);
    [Gt, dGt] = dop( startX, noBfs );
    Dt = dGt * Gt';
    %
    for k=1:length(halfRange)
        row = halfRange(k) * ones(length(range),1);
        %
        rows = [rows; row];
        cols = [cols; range];
        vals = [vals; Dt(halfRange(k),:)'];
    end;
    %
    % Compute the strip diagonal entries
    %
    noOnDiag = noPts - 2 * ls2;
    for k=1:noOnDiag
        localX = x(range+k);
        [Gt, dGt] = dop( localX, noBfs );
        tdGt = dGt(ls2,:);
        dt = tdGt * Gt';
        row = (k + ls2) * ones( length(range),1 );
        %
        rows = [rows; row];
        cols = [cols; range + k];
        vals = [vals; dt'];
    end;
    %
    % generate the bottom part of Si
    %
    endX = x(end-ls+1:end);
    [Gt, dGt] = dop( endX, noBfs );
    Dt = dGt * Gt';
    halfRange = (noPts-ls2+1:noPts)';
    range = (noPts-ls+1:noPts)';
    for k=1:length(halfRange)
        row = halfRange(k) * ones(length(range),1);
        %
        rows = [rows; row];
        cols = [cols; range];
        vals = [vals; Dt(k+ls2-1,:)'];
    end;
    %
    S = sparse(rows, cols, vals, noPts, noPts );
    %
    rS = sprank( S );
    if rS < noPts - 1
        warning(['The rank of S is ',int2str(rS),' while x has n = ',int2str(noPts),' points.']);
    end;
    %
    if ~genSparse
        S = full(S);
    end;
    %
end;
end

function [P,dP, recurrenceCoeffs] = dop( m, n )
%
% Function : Generates a set of discrete orthonormal polynomials, P, and
%   their derivatives, dP, either of size (m x n), or on the arbitrary
%   support vector x. It also returns the coefficients for the recurrence
%   relationships, these can be used to perform interpolation.
%
% Syntax :
%   [P,dP] = dop( m ) ;
%   [P,dP] = dop( m, n ) ;
%   [P,dP] = dop( x ) ;
%   [P,dP, rC] = dop( x, n ) ;
%
% Input :
%   m := number of evenly (unit) spaced points in support
%   x := arbitrary support vector
%   n := number of functions
%   
% Output :
%   P = [ p0, p1, ..., p(n-1) ] := Discete polynomials, pk are vectors.
%   dP = [ dp0, dp1, ..., dp(n-1) ] := Derivatives of the polynomials.
%   rC = a matrix, containing the alpha (col 1) and beta (col 2) 
%       coefficients for the three term recurrence relationship
%
% Cite this as :

% @article{DBLP:journals/tim/OLearyH12,
%  author    = {Paul O'Leary and
%               Matthew Harker},
%  title     = {A Framework for the Evaluation of Inclinometer Data in the
%               Measurement of Structures},
%  journal   = {IEEE T. Instrumentation and Measurement},
%  volume    = {61},
%  number    = {5},
%  year      = {2012},
%  pages     = {1237-1251},
% }
%
% @inproceedings{
% olearyHarker2008B,
%   Author = {O'Leary, Paul and Harker, Matthew},
%   Title = {An Algebraic Framework for Discrete Basis Functions in Computer Vision},
%   BookTitle = {IEEE Indian Conference on Computer Vision, Graphics and Image Processing},
%   Address= {Bhubaneswar, Dec},
%   Year = {2008} }
%
% Author : Matthew Harker
% Date : Nov. 29, 2011
% Version : 1.0
%--------------------------------------------------------------------------
% (c) 2011, Harker, O'Leary, University of Leoben, Leoben, Austria
% email: automation@unileoben.ac.at, url: automation.unileoben.ac.at
%--------------------------------------------------------------------------
% History:
%   Date:           Comment:
%   Nov. 29, 2011   Original Version 1.0
%--------------------------------------------------------------------------
%
[u,v] = size( m ) ;
%
if u == 1 && v == 1
    %
    x = (-1:2/(m-1):1)' ;
    %
elseif u ~= 1 && v == 1
    %
    x = m ;
    m = length(x) ;
    %
else
    %
    error('Support x should be an m x 1 vector') ;
    %
end
%
if nargin == 1
    n = m ;
end
%
%==============================
% Generate the Basis
%==============================
%
% Generate the first two polynomials :
p0 = ones(m,1)/sqrt(m) ;
meanX = mean( x );
p1 = x - meanX ;
np1 = norm( p1 ) ;
p1 = p1 / np1;
%
% Compute the derivatives of the degree-1 polynomial :
hm = sum( diff( x ) ) ; % Alternatively mean(...)
h = sum( diff( p1 ) ) ; % Alternatively mean(...), but 1/n cancels.
dp1 = (h/hm) * ones(m,1) ;
%
% Initialize the basis function matrices :
P = zeros(m,n) ;
P(:,1:2) = [ p0, p1 ] ;
%
dP = zeros(m,n) ;
dP(:,2) = dp1 ;
%
% Setup storage for the coefficients of the three term relationship
%
alphas = zeros(n,1);
alphas(1) = 1/sqrt(m);
alphas(2) = 1/np1;
%
betas = zeros(n,1);
betas(2) = meanX;
%
for k = 3:n
    %
    % Augment previous polynomial :
    pt = P(:,k-1) .* p1 ;
    %
    % 3-term recurrence :
    beta0 = (P(:,k-2)'*pt) ;
    pt = pt - P(:,k-2) *  beta0 ;
    betas(k) = beta0;
    %
    % Complete reorthogonalization :
    beta = P(:,1:k-1)' * pt ;
    pt = pt - P(:,1:k-1) * beta ;
    %
    % Apply coefficients to recurrence formulas : 
    alpha = 1/sqrt(pt'*pt) ;
    alphas(k) = alpha;
    P(:,k) = alpha * pt ;
    dP(:,k) = alpha * ( dP(:,k-1) .* p1 + P(:,k-1) .* dp1 - dP(:,k-2) * beta0 - dP(:,1:k-1)*beta  ) ;
    %
end;
%
recurrenceCoeffs = [alphas, betas];
%========
% END
%========
end

function even = isEven( numbers );
%
% Purpose : Checks is a number is even
%
% Use (syntax): even = isEven( numbers );
%
% Input Parameters :
%       numbers: a scalar or vector of scalars
%
% Return Parameters :
%       even: boolean or vector of boolean
%
% Description and algorithms:
%
% References : 
%
% Author :  Matthew Harker and Paul O'Leary
% Date :    17. January 2012
% Version : 1.0
%
% (c) 2013 Matthew Harker and Paul O'Leary, 
% Chair of Automation, University of Leoben, Leoben, Austria
% email: office@harkeroleary.org, 
% url: www.harkeroleary.org
%
% History:
%   Date:           Comment:
%

even = zeros( size( numbers ) );
%
rest = mod( numbers, 2 );
%
ind = find( rest == 0 );
%
if ~isempty( ind )
    even( ind ) = 1;
end;
end
