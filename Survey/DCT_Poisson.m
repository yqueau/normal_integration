function z = DCT_Poisson(p,q)
% An implementation of the use of DCT for solving the Poisson equation,
% (integration with Neumann boundary condition) 
% Code is based on the description in [1], Sec. 3.4
%
% [1] Normal Integration: a Survey - Queau et al., 2017
%
% Usage : 
% u=DCT_Poisson(p,q) 
% where p and q are MxN matrices, solves in the least square sense 
% \nabla u = [p,q] , assuming natural Neumann boundary condition
%
% \nabla u \cdot \eta = [p,q] \cdot \eta on boundaries
%
% Axis : O->y
%        |
%        x
%
% Fast solution is provided by Discrete Cosine Transform
%
% Implementation : Yvain Queau

% Divergence of (p,q) using central differences
px = 0.5*(p([2:end end],:)-p([1 1:end-1],:));
qy = 0.5*(q(:,[2:end end])-q(:,[1 1:end-1]));

% Div(p,q) 
f = px+qy;

% Right hand side of the boundary condition
b = zeros(size(p)); 
b(1,2:end-1) = -p(1,2:end-1);
b(end,2:end-1) = p(end,2:end-1);
b(2:end-1,1) = -q(2:end-1,1);
b(2:end-1,end) = q(2:end-1,end);
b(1,1) = (1/sqrt(2))*(-p(1,1)-q(1,1));
b(1,end) = (1/sqrt(2))*(-p(1,end)+q(1,end));
b(end,end) = (1/sqrt(2))*(p(end,end)+q(end,end));
b(end,1) = (1/sqrt(2))*(p(end,1)-q(end,1));

% Modification near the boundaries to enforce the non-homogeneous Neumann BC (Eq. 53 in [1])
f(1,2:end-1) = f(1,2:end-1)-b(1,2:end-1); 
f(end,2:end-1) = f(end,2:end-1)-b(end,2:end-1);
f(2:end-1,1) = f(2:end-1,1)-b(2:end-1,1); 
f(2:end-1,end) = f(2:end-1,end)-b(2:end-1,end);

% Modification near the corners (Eq. 54 in [1])
f(1,end) = f(1,end)-sqrt(2)*b(1,end);
f(end,end) = f(end,end)-sqrt(2)*b(end,end);
f(end,1) = f(end,1)-sqrt(2)*b(end,1);
f(1,1) = f(1,1)-sqrt(2)*b(1,1);

% Cosine transform of f
fcos=dct2(f);


% Cosine transform of z (Eq. 55 in [1])
[x,y] = meshgrid(0:size(p,2)-1,0:size(p,1)-1);
denom = 4*((sin(0.5*pi*x/size(p,2))).^2 + (sin(0.5*pi*y/size(p,1))).^2); 
z_bar_bar = -fcos./max(eps,denom);

% Inverse cosine transform :
z = idct2(z_bar_bar);
z=z-min(z(:)); % Z known up to a positive constant, so offset it to get from 0 to max

return

