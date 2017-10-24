clear
close all

addpath('Toolbox/');

% Tested methods
test_TV = 1; % Total variation
test_NC = 1; % Non-convex
test_AD = 1; % Anisotropic diffusion
test_MS = 1; % Mumford-Shah

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load a dataset containing:
% -- p (nrows x ncols) : gradient in the u- (bottom) direction
% -- q (nrows x ncols) : gradient in the v- (right) direction
% -- u (nrows x ncols) : ground truth depth map
% -- mask (nrows x ncols) : mask of the pixels on the vase (binary)
load Datasets/vase

% In this test we assume no mask is given, so discontinuities around the border should be recovered automatically
mask = ones(size(p));
indices_mask = find(mask>0);

% Add zero-mean, Gaussian noise
std_noise = 0.005*max(sqrt(p(indices_mask).^2+q(indices_mask).^2));
p = p+std_noise*randn(size(p));
q = q+std_noise*randn(size(q));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quadratic Integration
disp('Doing quadratic integration');

lambda = 1e-6*ones(size(p)); % Uniform field of weights (nrows x ncols)
z0 = zeros(size(p)); % Null depth prior (nrows x ncols)
solver = 'pcg'; % Solver ('pcg' means conjugate gradient, 'direct' means backslash i.e. sparse Cholesky) 
precond = 'CMG';  % Preconditioner for smooth integration ('none' means no preconditioning, 'ichol' means incomplete Cholesky, 'CMG' means conjugate combinatorial multigrid -- the latter is fastest, but it need being installed, see README)

t_1 = tic;
z_1 = smooth_integration(p,q,mask,lambda,z0,solver,precond);
t_1 = toc(t_1);

% Find the integration constant which minimizes RMSE
lambda_1 = -mean(z_1(indices_mask)-u(indices_mask));
z_1 = z_1+lambda_1;
% Calculate RMSE
RMSE_1 = sqrt(mean((z_1(indices_mask)-u(indices_mask)).^2));
% Display evaluation results in terminal
disp('=============================');
disp('Quadratic integration:');
disp(sprintf('CPU: %.4f',t_1)); 
disp(sprintf('RMSE:   %.2f',RMSE_1));
disp(' ');

if(test_TV)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% TV Integration
	disp('Doing TV integration');

	zinit = z_1; % least-squares initialization
	alpha = 0.1; % Descent stepsize (influences speed)
	tol = 1e-5; % Stopping criterion
	maxit = 1000; % Stopping criterion

	t_2 = tic;
	z_2 = tv_integration(p,q,mask,lambda,z0,alpha,maxit,tol,zinit);
	t_2 = toc(t_2);

	% Find the integration constant which minimizes RMSE
	lambda_2 = -mean(z_2(indices_mask)-u(indices_mask));
	z_2 = z_2+lambda_2;
	% Calculate RMSE
	RMSE_2 = sqrt(mean((z_2(indices_mask)-u(indices_mask)).^2));
	% Display evaluation results in terminal
	disp('=============================');
	disp('TV integration:');
	disp(sprintf('CPU: %.4f',t_2)); 
	disp(sprintf('RMSE:   %.2f',RMSE_2));
	disp(' ');
end

if(test_NC)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Nonconvex Integration
	disp('Doing nonconvex integration');

	zinit = z_1; % least-squares initialization
	gamma = 0.5; % Nonconvex estimator parameter (to be tuned: e.g. 0.5 for phi1, 1 for phi2 in our tests)
	beta = 0.8; % Lischitz reduction constant (must be in (0,1), see iPiano paper, 0.8 seems to always work)
	maxit = 1000; % Stopping criterion
	tol = 1e-5; % Stopping criterion
	
	t_3 = tic;
	z_3 = phi1_integration(p,q,mask,lambda,z0,beta,gamma,maxit,tol,zinit,u); % Phi_1 estimator 
	% z_3 = phi2_integration(p,q,mask,lambda,z0,beta,gamma,maxit,tol,zinit,u);%  Phi_2 estimator 
	t_3 = toc(t_3);

	% Find the integration constant which minimizes RMSE
	lambda_3 = -mean(z_3(indices_mask)-u(indices_mask));
	z_3 = z_3+lambda_3;
	% Calculate RMSE
	RMSE_3 = sqrt(mean((z_3(indices_mask)-u(indices_mask)).^2));
	% Display evaluation results in terminal
	disp('=============================');
	disp('Nonconvex integration:');
	disp(sprintf('CPU: %.4f',t_3)); 
	disp(sprintf('RMSE:   %.2f',RMSE_3));
	disp(' ');
end

if(test_AD)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Anisotropic diffusion ontegration
	disp('Doing anis diff integration');

	zinit = z_1; % least-squares initialization
	mu = 0.2;	% anis diff (to be tuned)
	nu = 10; % anis diff param (10 should work)
	maxit = 20; % Stopping criterion
	tol = 1e-5; % Stopping criterion
	
	t_4 = tic;
	z_4 = anisotropic_diffusion_integration(p,q,mask,lambda,z0,mu,nu,maxit,tol,zinit);
	t_4 = toc(t_4);

	% Find the integration constant which minimizes RMSE
	lambda_4 = -mean(z_4(indices_mask)-u(indices_mask));
	z_4 = z_4+lambda_4;
	% Calculate RMSE
	RMSE_4 = sqrt(mean((z_4(indices_mask)-u(indices_mask)).^2));
	% Display evaluation results in terminal
	disp('=============================');
	disp('Anis diff integration:');
	disp(sprintf('CPU: %.4f',t_4)); 
	disp(sprintf('RMSE:   %.2f',RMSE_4));
	disp(' ');
end

if(test_MS)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Mumford-Shah ontegration
	disp('Doing Mumford-Shah integration');

	zinit = z_1; % least-squares initialization
	mu = 45; % Regularization weight for discontinuity set
	epsilon = 0.01; % Should be close to 0
	tol = 1e-5; % Stopping criterion
	maxit = 1000; % Stopping criterion 
	
	t_5 = tic;
	z_5 = mumford_shah_integration(p,q,mask,lambda,z0,mu,epsilon,maxit,tol,zinit);
	t_5 = toc(t_5);

	% Find the integration constant which minimizes RMSE
	lambda_5 = -mean(z_5(indices_mask)-u(indices_mask));
	z_5 = z_5+lambda_5;
	% Calculate RMSE
	RMSE_5 = sqrt(mean((z_5(indices_mask)-u(indices_mask)).^2));
	% Display evaluation results in terminal
	disp('=============================');
	disp('Mumford-Shah integration:');
	disp(sprintf('CPU: %.4f',t_5)); 
	disp(sprintf('RMSE:   %.2f',RMSE_5));
	disp(' ');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summarize results
disp(' ');
disp(' ');
disp(' ');
disp('=============================');
disp('=============================');
disp('Summary of the evaluation:');
disp('=============================');
disp('=============================');
disp('Quadratic integration:');
disp(sprintf('CPU: %.4f',t_1)); 
disp(sprintf('RMSE:   %.2f',RMSE_1));
disp(' ');
if(test_TV)
	% Display evaluation results in terminal
	disp('=============================');
	disp('TV integration:');
	disp(sprintf('CPU: %.4f',t_2)); 
	disp(sprintf('RMSE:   %.2f',RMSE_2));
	disp(' ');
end
if(test_NC)
	% Display evaluation results in terminal
	disp('=============================');
	disp('Nonconvex integration:');
	disp(sprintf('CPU: %.4f',t_3)); 
	disp(sprintf('RMSE:   %.2f',RMSE_3));
	disp(' ');
end
if(test_AD)
	% Display evaluation results in terminal
	disp('=============================');
	disp('Anis diff integration:');
	disp(sprintf('CPU: %.4f',t_4)); 
	disp(sprintf('RMSE:   %.2f',RMSE_4));
	disp(' ');
end
if(test_MS)
	% Display evaluation results in terminal
	disp('=============================');
	disp('Mumford-Shah integration:');
	disp(sprintf('CPU: %.4f',t_5)); 
	disp(sprintf('RMSE:   %.2f',RMSE_5));
	disp(' ');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display a few things

figure('units','normalized','outerposition',[0 0 1 1])

% Input data: p, q and mask
subplot(4,4,1)
imagesc(p);
axis image
axis off
title('$$p$$','Interpreter','Latex','Fontsize',14)
subplot(4,4,2)
imagesc(q);
axis image
axis off
title('$$q$$','Interpreter','Latex','Fontsize',14)
subplot(4,4,3)
imagesc(mask);
axis image
axis off
colormap gray
title('$$\Omega$$','Interpreter','Latex','Fontsize',14)
subplot(4,4,4)
surfl(u,[-135 30]);
view(-35,20)
axis ij;
shading flat;
colormap gray;
axis equal;
grid off
axis([1 size(p,1) 1 size(p,2) min(u(:)) max(u(:))]);
axis off
title('Ground truth depth','Interpreter','Latex','Fontsize',14)


subplot(4,4,5)
surfl(z_1,[-135 30]);
view(-35,20)
axis ij;
shading flat;
colormap gray;
axis equal;
grid off
axis([1 size(p,1) 1 size(p,2) min(u(:)) max(u(:))]);
axis off
title('Quadratic integration','Interpreter','Latex','Fontsize',14)

error_map_1 = abs(u-z_1);
error_map_1(mask==0) = NaN;

subplot(4,4,6)
imagesc(error_map_1,[0 10]);
axis image
axis off
colormap gray
title('Absolute error (quadratic integration)','Interpreter','Latex','Fontsize',14)

if(test_TV)
	subplot(4,4,7)
	surfl(z_2,[-135 30]);
	view(-35,20)
	axis ij;
	shading flat;
	colormap gray;
	axis equal;
	grid off
	axis([1 size(p,1) 1 size(p,2) min(u(:)) max(u(:))]);
	axis off
	title('TV integration','Interpreter','Latex','Fontsize',14)

	error_map_2 = abs(u-z_2);
	error_map_2(mask==0) = NaN;

	subplot(4,4,8)
	imagesc(error_map_2,[0 10]);
	axis image
	axis off
	colormap gray
	title('Absolute error (TV integration)','Interpreter','Latex','Fontsize',14)
end


if(test_NC)
	subplot(4,4,9)
	surfl(z_3,[-135 30]);
	view(-35,20)
	axis ij;
	shading flat;
	colormap gray;
	axis equal;
	grid off
	axis([1 size(p,1) 1 size(p,2) min(u(:)) max(u(:))]);
	axis off
	title('Nonconvex integration','Interpreter','Latex','Fontsize',14)

	error_map_3 = abs(u-z_3);
	error_map_3(mask==0) = NaN;

	subplot(4,4,10)
	imagesc(error_map_3,[0 10]);
	axis image
	axis off
	colormap gray
	title('Absolute error (nonconvex integration)','Interpreter','Latex','Fontsize',14)
end

if(test_AD)
	subplot(4,4,11)
	surfl(z_4,[-135 30]);
	view(-35,20)
	axis ij;
	shading flat;
	colormap gray;
	axis equal;
	grid off
	axis off
	axis([1 size(p,1) 1 size(p,2) min(u(:)) max(u(:))]);
	title('Anis diff integration','Interpreter','Latex','Fontsize',14)

	error_map_4 = abs(u-z_4);
	error_map_4(mask==0) = NaN;

	subplot(4,4,12)
	imagesc(error_map_4,[0 10]);
	axis image
	axis off
	colormap gray
	title('Absolute error (anis diff integration)','Interpreter','Latex','Fontsize',14)
end


if(test_MS)
	subplot(4,4,13)
	surfl(z_5,[-135 30]);
	view(-35,20)
	axis ij;
	shading flat;
	colormap gray;
	axis equal;
	grid off
	axis([1 size(p,1) 1 size(p,2) min(u(:)) max(u(:))]);
	axis off
	title('Mumford-Shah integration','Interpreter','Latex','Fontsize',14)

	error_map_5 = abs(u-z_5);
	error_map_5(mask==0) = NaN;

	subplot(4,4,14)
	imagesc(error_map_5,[0 10]);
	axis image
	axis off
	colormap gray
	title('Absolute error (Mumford-Shah integration)','Interpreter','Latex','Fontsize',14)
end
