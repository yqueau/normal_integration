clear all
close all

% Noise standard deviations which are being tested
tab_std = 0:0.025:0.25;

% Structure to store the RMSE of all three methods
tab_RMSE = zeros(length(tab_std),1);
tab_RMSE_HOL = zeros(length(tab_std),1);
tab_RMSE_SCS = zeros(length(tab_std),1);
cpt = 1;

% For each size, do nb_trials tests and mean the result to reduce variability of measurements
nb_trials = 50;
tab_local = zeros(nb_trials,1);
tab_local_HOL = zeros(nb_trials,1);
tab_local_SCS = zeros(nb_trials,1);

% Create a m x n dataset
n = 64;
m = n;
[u,q_p,p_p,y,x] = g2sTestSurf(n,m,'even',0);  
% Data is sampled between -1 and 10 : rescale it
facteur = 11/n;
u = u/facteur;
% Mask (here, full rectangular domain)
mask = ones(size(p_p));
imask = find(mask>0);

% Run the evaluation
for l = tab_std

	for trial = 1:nb_trials

		% Add noise
		std_noise = l*max(sqrt(p_p(imask).^2+q_p(imask).^2));
		p = p_p+std_noise*randn(size(p_p));
		q = q_p+std_noise*randn(size(p_p));

		% Solve by ours
		z = smooth_integration(p,q,mask,[],[],'pcg','CMG');
		% Find integration constant which minimizes RMSE
		lambda = -mean(z(imask)-u(imask));
		z = z+lambda;
		% Evaluate
		RMSE = sqrt(mean((z(imask)-u(imask)).^2));
		tab_local(trial) = RMSE;

		% Solve by SCS
		z = simchony(p,q);
		% Find integration constant which minimizes RMSE
		lambda = -mean(z(imask)-u(imask));
		z = z+lambda;
		% Evaluate
		RMSE = sqrt(mean((z(imask)-u(imask)).^2));
		tab_local_SCS(trial) = RMSE;

		% Solve by HOL
		z = harker_oleary(p,q);
		% Find integration constant which minimizes RMSE
		lambda = -mean(z(imask)-u(imask));
		z = z+lambda;
		% Evaluate
		RMSE = sqrt(mean((z(imask)-u(imask)).^2));
		tab_local_HOL(trial) = RMSE;

	end

	RMSE = mean(tab_local);
	RMSE_SCS = mean(tab_local_SCS);
	RMSE_HOL = mean(tab_local_HOL);

	% Display
	disp(sprintf('sigma = %.2f - RMSE = %.4f',l,RMSE));
	disp(sprintf('sigma = %.2f - RMSE - SCS = %.4f',l,RMSE_SCS));
	disp(sprintf('sigma = %.2f - RMSE - HOL = %.4f',l,RMSE_HOL));

	tab_RMSE(cpt) = RMSE;
	tab_RMSE_SCS(cpt) = RMSE_SCS;
	tab_RMSE_HOL(cpt) = RMSE_HOL;

	cpt = cpt+1;
end

figure(1)
plot(tab_std,tab_RMSE_SCS,'b-.','Linewidth',4)
hold on
plot(tab_std,tab_RMSE_HOL,'r--','Linewidth',4)
plot(tab_std,tab_RMSE,'k-','Linewidth',4)
set(gca,'FontSize',18)
xlabel('$$\sigma$$','Interpreter','Latex','Fontsize',24);
ylabel('RMSE ($px$)','Interpreter','Latex','Fontsize',24);
hl = legend('Simchony et al.','Harker and O''Leary','Proposed');
set(hl,'Interpreter','Latex')
