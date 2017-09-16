clear
close all

addpath('Toolbox/');

% Tested methods
test_FFT = 1; % FFT integrator (periodic BC)
test_DST = 1; % DST integrator (Dirichlet BC)
test_DCT = 1; % DCT integrator (natural BC)
test_HB = 1; % Modified Horn and Brook's scheme

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load a dataset containing:
% -- p (nrows x ncols) : gradient in the u- (bottom) direction
% -- q (nrows x ncols) : gradient in the v- (right) direction
% -- u (nrows x ncols) : ground truth depth map
% -- mask (nrows x ncols) : mask of the pixels on the vase (binary)
load Datasets/vase

% To emphasize the problem of boundaries, we crop the domain t
p = p(83:310,90:180);
q = q(83:310,90:180);
u = u(83:310,90:180);
mask = mask(83:310,90:180);
%~ p = p(83:260,150:220);
%~ q = q(83:260,150:220);
%~ u = u(83:260,150:220);
%~ mask = mask(83:260,150:220);
indices_mask = find(mask>0);

% Add zero-mean, Gaussian noise
std_noise = 0.005*max(sqrt(p(indices_mask).^2+q(indices_mask).^2));
p(indices_mask) = p(indices_mask)+std_noise*randn(size((indices_mask)));
q(indices_mask) = q(indices_mask)+std_noise*randn(size((indices_mask)));

if(test_FFT)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% FFT Integration
	disp('Doing FFT integration');


	t_1 = tic;
	z_1 = FFT_Poisson(p,q);
	t_1 = toc(t_1);

	% Find the integration constant which minimizes RMSE
	lambda_1 = -mean(z_1(indices_mask)-u(indices_mask));
	z_1 = z_1+lambda_1;
	% Calculate RMSE
	RMSE_1 = sqrt(mean((z_1(indices_mask)-u(indices_mask)).^2));
	% Display evaluation results in terminal
	disp('=============================');
	disp('FFT integration:');
	disp(sprintf('CPU: %.4f',t_1)); 
	disp(sprintf('RMSE:   %.2f',RMSE_1));
	disp(' ');
end

if(test_DST)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DST Integration
	disp('Doing DST integration');

	u_b = zeros(size(p)); % Homogeneous Dirichlet BC
	% u_b(1,:) = 100; % For a more funny boundary, uncomment this ;)
	t_2 = tic;
	z_2 = DST_Poisson(p,q,u_b);
	t_2 = toc(t_2);

	% Find the integration constant which minimizes RMSE
	lambda_2 = -mean(z_2(indices_mask)-u(indices_mask));
	z_2 = z_2+lambda_2;
	% Calculate RMSE
	RMSE_2 = sqrt(mean((z_2(indices_mask)-u(indices_mask)).^2));
	% Display evaluation results in terminal
	disp('=============================');
	disp('DST integration:');
	disp(sprintf('CPU: %.4f',t_2)); 
	disp(sprintf('RMSE:   %.2f',RMSE_2));
	disp(' ');
end

if(test_DCT)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DCT Integration
	disp('Doing DCT integration');

	t_3 = tic;
	z_3 = DCT_Poisson(p,q); % Natural Neumann BC
	t_3 = toc(t_3);

	% Find the integration constant which minimizes RMSE
	lambda_3 = -mean(z_3(indices_mask)-u(indices_mask));
	z_3 = z_3+lambda_3;
	% Calculate RMSE
	RMSE_3 = sqrt(mean((z_3(indices_mask)-u(indices_mask)).^2));
	% Display evaluation results in terminal
	disp('=============================');
	disp('DCT integration:');
	disp(sprintf('CPU: %.4f',t_3)); 
	disp(sprintf('RMSE:   %.2f',RMSE_3));
	disp(' ');
end

if(test_HB)
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Anisotropic diffusion ontegration
	disp('Doing Horn Brooks integration');

	zinit = z_1; % least-squares initialization
	maxit = 50000; % Stopping criterion
	tol = 1e-6; % Stopping criterion
	trace = 0; % To display or not the current estimate
	
	t_4 = tic;
	z_4 = horn_brooks(p,q,mask,maxit,tol,trace);
	t_4 = toc(t_4);

	% Find the integration constant which minimizes RMSE
	lambda_4 = -mean(z_4(indices_mask)-u(indices_mask));
	z_4 = z_4+lambda_4;
	% Calculate RMSE
	RMSE_4 = sqrt(mean((z_4(indices_mask)-u(indices_mask)).^2));
	% Display evaluation results in terminal
	disp('=============================');
	disp('Horn-Brook integration:');
	disp(sprintf('CPU: %.4f',t_4)); 
	disp(sprintf('RMSE:   %.2f',RMSE_4));
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

if(test_FFT)
	disp('=============================');
	disp('FFT integration:');
	disp(sprintf('CPU: %.4f',t_1)); 
	disp(sprintf('RMSE:   %.2f',RMSE_1));
	disp(' ');
end
if(test_DST)
	% Display evaluation results in terminal
	disp('=============================');
	disp('DST integration:');
	disp(sprintf('CPU: %.4f',t_2)); 
	disp(sprintf('RMSE:   %.2f',RMSE_2));
	disp(' ');
end
if(test_DCT)
	% Display evaluation results in terminal
	disp('=============================');
	disp('DCT integration:');
	disp(sprintf('CPU: %.4f',t_3)); 
	disp(sprintf('RMSE:   %.2f',RMSE_3));
	disp(' ');
end
if(test_HB)
	% Display evaluation results in terminal
	disp('=============================');
	disp('Horn Brooks integration:');
	disp(sprintf('CPU: %.4f',t_4)); 
	disp(sprintf('RMSE:   %.2f',RMSE_4));
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
view(-60,20)
axis ij;
shading flat;
colormap gray;
axis equal;
grid off
axis([1 size(p,2) 1 size(p,1) min(u(:)) max(u(:))]);
axis off
title('Ground truth depth','Interpreter','Latex','Fontsize',14)


if(test_FFT)
	subplot(4,4,5)
	surfl(z_1,[-135 30]);
	view(-60,20)
	axis ij;
	shading flat;
	colormap gray;
	axis equal;
	grid off
	axis([1 size(p,2) 1 size(p,1) min(u(:)) max(u(:))]);
	axis off
	title('FFT integration','Interpreter','Latex','Fontsize',14)

	error_map_1 = abs(u-z_1);
	error_map_1(mask==0) = NaN;

	subplot(4,4,6)
	imagesc(error_map_1,[0 10]);
	axis image
	axis off
	colormap gray
	title('Absolute error (FFT integration)','Interpreter','Latex','Fontsize',14)
end

if(test_DST)
	subplot(4,4,7)
	surfl(z_2,[-135 30]);
	view(-60,20)
	axis ij;
	shading flat;
	colormap gray;
	axis equal;
	grid off
	axis([1 size(p,2) 1 size(p,1) min(u(:)) max(u(:))]);
	axis off
	title('DST integration','Interpreter','Latex','Fontsize',14)

	error_map_2 = abs(u-z_2);
	error_map_2(mask==0) = NaN;

	subplot(4,4,8)
	imagesc(error_map_2,[0 10]);
	axis image
	axis off
	colormap gray
	title('Absolute error (DST integration)','Interpreter','Latex','Fontsize',14)
end


if(test_DCT)
	subplot(4,4,9)
	surfl(z_3,[-135 30]);
	view(-60,20)
	axis ij;
	shading flat;
	colormap gray;
	axis equal;
	grid off
	axis([1 size(p,2) 1 size(p,1) min(u(:)) max(u(:))]);
	axis off
	title('DCT integration','Interpreter','Latex','Fontsize',14)

	error_map_3 = abs(u-z_3);
	error_map_3(mask==0) = NaN;

	subplot(4,4,10)
	imagesc(error_map_3,[0 10]);
	axis image
	axis off
	colormap gray
	title('Absolute error (DCT integration)','Interpreter','Latex','Fontsize',14)
end

if(test_HB)
	subplot(4,4,11)
	surfl(z_4,[-135 30]);
	view(-60,20)
	axis ij;
	shading flat;
	colormap gray;
	axis equal;
	grid off
	axis off
	axis([1 size(p,2) 1 size(p,1) min(u(:)) max(u(:))]);
	title('Horn Brooks integration','Interpreter','Latex','Fontsize',14)

	error_map_4 = abs(u-z_4);
	error_map_4(mask==0) = NaN;

	subplot(4,4,12)
	imagesc(error_map_4,[0 10]);
	axis image
	axis off
	colormap gray
	title('Absolute error (Horn Brooks integration)','Interpreter','Latex','Fontsize',14)
end

