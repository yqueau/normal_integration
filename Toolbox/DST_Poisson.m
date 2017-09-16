function z = DST_Poisson(p,q,ub)
% An implementation of the use of DST for solving the Poisson equatio,
% (integration with Dirichlet boundary condition) 
% Code is based on the description in [1], Sec. 3.4
%
% [1] Normal Integration: a Survey - Queau et al., 2017
%
% Usage : 
% u=DST_Poisson(p,q) 
% where p and q are MxN matrices, solves in the least square sense 
% \nabla u = [p,q] , assuming homogeneous Dirichlet boundary 
% condition u = 0 on the boundary
%
% u=DST_Poisson(p,q,ub) 
% where p,q and ub are NxM matrix, such as ub(1,:) contains the values 
% on the first line, ub(end,:) contains the values on the last line etc.
% Apart from the boundary ub can be anything
%
% Example (weird)
% p = zeros(100,100); 
% q = zeros(100,100);
% ub = zeros(size(p));
% ub(1,:)=1:1OO; 
% u = DST_Poisson(p,q,ub);
% surfl(u)
%
% This performs the least square solution to \nabla u = [p,q], i.e. :
% min \int_\Omega \| \nablua U - [p,q] \|^2
% where \Omega is square and the Dirichlet boundary condition 
% u = ub on the boundary of \Omega. 
%
% Axis : O->y
%        |
%        x
%
% Fast solution is provided by Discrete Sine Transform
%
% Implementation : Yvain Queau


if(nargin<3)
	ub=zeros(size(p));
end

% Divergence of (p,q) using central differences
px = 0.5*(p([2:end end],:)-p([1 1:end-1],:));
qy = 0.5*(q(:,[2:end end])-q(:,[1 1:end-1]));
f = px + qy;

% Modification near the boundaries (Eq. 46 in [1])
f(2,3:end-2) = f(2,3:end-2) - ub(1,3:end-2); 
f(end-1,3:end-2) = f(end-1,3:end-2) - ub(end,3:end-2); 
f(3:end-2,2) = f(3:end-2,2) - ub(3:end-2,1); 
f(3:end-2,end-1) = f(3:end-2,end-1) - ub(3:end-2,end);

% Modification near the corners (Eq. 47 in [1])
f(2,2) = f(2,2) - ub(2,1) - ub(1,2); 
f(2,end-1) = f(2,end-1) - ub(2,end) - ub(1,end-1); 
f(end-1,end-1) = f(end-1,end-1) - ub(end-1,end) - ub(end,end-1); 
f(end-1,2) = f(end-1,2) - ub(end-1,1) - ub(end,2);

% Sine transform of f
fsin=dst2(f(2:end-1,2:end-1));

% Denominator
[x,y] = meshgrid(0:size(p,2)-1,0:size(p,1)-1);
denom = (sin(0.5*pi*x/size(p,2))).^2 + (sin(0.5*pi*y/size(p,1))).^2; 
z_bar = -0.25*fsin./denom(2:end-1,2:end-1);

% Inverse Sine transform :
z=ub;
z(2:end-1,2:end-1) = idst2(z_bar);

return


function y = dst2(x)
y = dst(dst(x)')';
return

function Y=idst2(X);
Z=idst(X');
Y=idst(Z');
return


