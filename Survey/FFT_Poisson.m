function z = FFT_Poisson(p,q,ub)
% An implementation of the use of FFT for solving the Poisson equation,
% (integration with periodic boundary condition) 
% Code is based on the description in [1], Sec. 3.3
%
% [1] Normal Integration: a Survey - Queau et al., 2017
%
% Usage : 
% u=FFT_Poisson(p,q) 
% where p and q are MxN matrices, solves in the least square sense 
% \nabla u = [p,q] , assuming periodic boundary condition 
%
% This performs the least square solution to \nabla u = [p,q], i.e. :
% min \int_\Omega \| \nablua U - [p,q] \|^2
% where \Omega is square and periodic boundary condition is enforced
%
% Axis : O->y
%        |
%        x
%
% Fast solution is provided by Fast Fourier Transform
%
% Implementation : Yvain Queau


% Fourier transforms of p and q
p_hat = fft2(p);
q_hat = fft2(q);

% Fourier transform of z (Eq. 42 in [1])
[y,x] = meshgrid(0:size(p,2)-1,0:size(p,1)-1);
numerator = sin(2*pi*x/size(p,1)).*p_hat+sin(2*pi*y/size(p,2)).*q_hat;
denominator = max(eps,4*j*((sin(pi*x/size(p,1))).^2+(sin(pi*y/size(p,2))).^2));

% Inverse Fourier transform
z = real(ifft2(numerator./denominator));
z=z-min(z(:)); % Z known up to a positive constant, so offset it to get from 0 to max

return



