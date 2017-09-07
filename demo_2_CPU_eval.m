clear all
close all

% Size of the domains which are being tested
tab_size = [128 256 512 1024 2048 4096];

% Structure to store the CPU times of all three methods
tab_CPU = zeros(length(tab_size),1);
tab_CPU_HOL = zeros(length(tab_size),1);
tab_CPU_SCS = zeros(length(tab_size),1);

% For each size, do nb_trials tests and mean the result to reduce variability of measurements
nb_trials = 5;
tab_local = zeros(nb_trials,1);
tab_local_HOL = zeros(nb_trials,1);
tab_local_SCS = zeros(nb_trials,1);

% Run the test
cpt = 1;
for n = tab_size

	for trial = 1:nb_trials
		% Make surface
		m = n;
		[u,q,p] = g2sTestSurf(n,m,'even',0);  
		facteur = 11/n;
		u = u/facteur;
		mask = ones(size(p));
		imask = find(mask>0);

		% Solve by ours
		tic
		z = smooth_integration(p,q,mask,[],[],'pcg','CMG');
		tab_local(trial) = toc;

		% Solve by SCS
		tic
		z = simchony(p,q);
		tab_local_SCS(trial) = toc;

		% Solve by HOL
		tic;
		z = harker_oleary(p,q);
		tab_local_HOL(trial) = toc;

	end

	CPU = mean(tab_local);
	CPU_SCS = mean(tab_local_SCS);
	CPU_HOL = mean(tab_local_HOL);

	% Display
	disp(sprintf('n = %d - CPU = %.4f',n,CPU));
	disp(sprintf('n = %d - CPU - SCS = %.4f',n,CPU_SCS));
	disp(sprintf('n = %d - CPU - HOL = %.4f',n,CPU_HOL));

	tab_CPU(cpt) = CPU;
	tab_CPU_SCS(cpt) = CPU_SCS;
	tab_CPU_HOL(cpt) = CPU_HOL;

	cpt = cpt+1;
end

% Show the CPU times
figure(1)
loglog(tab_size.^2,tab_CPU_SCS,'b-.','Linewidth',4)
hold on
loglog(tab_size.^2,tab_CPU_HOL,'r--','Linewidth',4)
loglog(tab_size.^2,tab_CPU,'k-','Linewidth',4)
hold off
set(gca,'FontSize',18)
xlabel('$$|\Omega|$$','Interpreter','Latex','Fontsize',24);
ylabel('CPU $$(s)$$','Interpreter','Latex','Fontsize',24);
hl = legend('Simchony et al.','Harker and O''Leary','Proposed');
set(hl,'Interpreter','Latex')
ax = axis;
axis([0,tab_size(end)^2,1e-3,1e3]);
set(gca,'XTick',tab_size.^2)
set(gca,'XTickLabel',{'128x128';'256x256';'512x512';'1024x1024';'2048x2048';'4096x4096'})
rotateXLabels( gca(), 45 );

