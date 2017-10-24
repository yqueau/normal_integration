clear
close all

addpath('Toolbox/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load a dataset containing:
% -- p (nrows x ncols) : gradient in the u- (bottom) direction
% -- q (nrows x ncols) : gradient in the v- (right) direction
% -- u (nrows x ncols) : ground truth depth map
% -- mask (nrows x ncols) : mask of the pixels on the vase (binary)
load Datasets/vase
indices_mask = find(mask>0); % Indices of the pixel inside the mask

% Add zero-mean, Gaussian noise inside the mask
std_noise = 0.02*max(sqrt(p(indices_mask).^2+q(indices_mask).^2));
p = p+std_noise*randn(size(p));
q = q+std_noise*randn(size(q));

% Fill the gradient with 0 to test rectangular integration
p(mask==0) = 0;
q(mask==0) = 0;
% Remove the ground truth depth values outside the mask
u(mask==0) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set optimization parameters
lambda = 1e-6*ones(size(p)); % Uniform field of weights (nrows x ncols)
z0 = zeros(size(p)); % Null depth prior (nrows x ncols)
solver = 'pcg'; % Solver ('pcg' means conjugate gradient, 'direct' means backslash i.e. sparse Cholesky) 
precond = 'CMG'; % Preconditioner ('none' means no preconditioning, 'ichol' means incomplete Cholesky, 'CMG' means conjugate combinatorial multigrid -- the latter is fastest, but it need being installed, see README)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate without mask (full rectangular domain)
Omega_1 = ones(size(p));
t_1 = tic;
z_1 = smooth_integration(p,q,Omega_1,lambda,z0,solver,precond);
t_1 = toc(t_1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integrate with mask (reduced non-rectangular domain)
Omega_2 = mask;
t_2 = tic;
z_2 = smooth_integration(p,q,Omega_2,lambda,z0,solver,precond);
t_2 = toc(t_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate both integration methods over the mask

% Find the constant of integration which minimizes RMSE wrt ground truth
lambda_1 = -mean(z_1(indices_mask)-u(indices_mask));
z_1 = z_1+lambda_1;

lambda_2 = -mean(z_2(indices_mask)-u(indices_mask));
z_2 = z_2+lambda_2;

% Calculate RMSEs 
RMSE_1 = sqrt(mean((z_1(indices_mask)-u(indices_mask)).^2));
RMSE_2 = sqrt(mean((z_2(indices_mask)-u(indices_mask)).^2));

% Display evaluation results in terminal
disp('=============================');
disp('Integration without mask:');
disp(sprintf('CPU: %.4f',t_1)); 
disp(sprintf('RMSE over mask:   %.2f',RMSE_1));
disp('');
disp('================================');
disp('Integration with mask:');
disp(sprintf('CPU: %.4f',t_2)); 
disp(sprintf('RMSE over mask:   %.2f',RMSE_2));
disp(''); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display a few things

figure('units','normalized','outerposition',[0 0 1 1])

% Input data: p, q and mask
subplot(3,3,1)
imagesc(p);
axis image
axis off
title('$$p$$','Interpreter','Latex','Fontsize',14)
subplot(3,3,2)
imagesc(q);
axis image
axis off
title('$$q$$','Interpreter','Latex','Fontsize',14)
subplot(3,3,3)
imagesc(mask);
axis image
axis off
colormap gray
title('$$\Omega$$','Interpreter','Latex','Fontsize',14)

subplot(3,3,4)
surfl(u,[-135 30]);
view(-35,20)
axis ij;
shading flat;
colormap gray;
axis equal;
grid off
axis off
title('Ground truth depth','Interpreter','Latex','Fontsize',14)


subplot(3,3,5)
surfl(z_1,[-135 30]);
view(-35,20)
axis ij;
shading flat;
colormap gray;
axis equal;
grid off
axis off
title('Reconstruction without mask','Interpreter','Latex','Fontsize',14)

subplot(3,3,6)
surfl(z_2,[-135 30]);
view(-35,20)
axis ij;
shading flat;
colormap gray;
axis equal;
grid off
axis off
title('Reconstruction with mask','Interpreter','Latex','Fontsize',14)

error_map_1 = abs(u-z_1);
error_map_1(mask==0) = NaN;

error_map_2 = abs(u-z_2);
error_map_2(mask==0) = NaN;

subplot(3,3,8)
imagesc(error_map_1,[0 5]);
axis image
axis off
colormap gray
title('Absolute error without mask','Interpreter','Latex','Fontsize',14)

subplot(3,3,9)
imagesc(error_map_2,[0 5]);
axis image
axis off
colormap gray
title('Absolute error with mask','Interpreter','Latex','Fontsize',14)
