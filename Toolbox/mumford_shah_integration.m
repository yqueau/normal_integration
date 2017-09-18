function [z,wup,wum,wvp,wvm,tab_nrj,tab_rmse] = mumford_shah_integration(p,q,mask,lambda,z0,mu,epsilon,maxit,tol,zinit,gt)

	% Check arguments
	if(nargin < 2)
		disp('Error: Not enough arguments');
		return;
	end
	
	[nrows,ncols] = size(p);
	
	% Set default values for missing arguments
	if (~exist('mask','var')|isempty(mask)) mask=ones(nrows,ncols); end;
	if (~exist('lambda','var')|isempty(lambda)) lambda = 1e-9*mask; end;
	if (~exist('z0','var')|isempty(z0)) z0 = zeros(nrows,ncols); end;
	if (~exist('mu','var')|isempty(mu)) mu = 1e-3; end;
	if (~exist('epsilon','var')|isempty(epsilon)) epsilon = 1e-3; end;
	if (~exist('maxit','var')|isempty(maxit)) maxit = 100; end;
	if (~exist('tol','var')|isempty(tol)) tol = 1e-4; end;
	if (~exist('zinit','var')|isempty(zinit)) zinit = z0; end;
	if (~exist('gt','var')|isempty(gt)) gt = z0; end;
	

	% If lambda is a scalar, make it a matrix
	if(size(lambda,1)==1)
		lambda = lambda*mask;
	end

	% Make finite differences operators	
	[Dup,Dum,Dvp,Dvm,imask,imaskup,imaskum,imaskvp,imaskvm] = make_gradient(mask);
	npix = length(imask);

	% Some stuff used later
	L = 0.5*(Dup'*Dup+Dum'*Dum+Dvp'*Dvp+Dvm'*Dvm); % Laplacian
	Lambda_two = spdiags(lambda(imask),0,npix,npix); % Regularization
	DuptDup = Dup'*Dup;
	DumtDum = Dum'*Dum;
	DvptDvp = Dvp'*Dvp;
	DvmtDvm = Dvm'*Dvm;
	bw = 0.5*((1/(4*epsilon))*ones(npix,1));
	if(nargout>5)
		tab_nrj = zeros(maxit+1,1);
		cpt = 1;
	end
	if(nargout>6)
		tab_rmse = zeros(maxit+1,1);
		cpt_rmse = 1;
	end
	

	% Initialization
	z = zinit;
	
	wup = zeros(npix,1);wup(imaskup) = 1;
	wum = zeros(npix,1);wum(imaskum) = 1;
	wvp = zeros(npix,1);wvp(imaskvp) = 1;
	wvm = zeros(npix,1);wvm(imaskvm) = 1;
	Wup2_mat = spdiags(wup.^2,0,npix,npix);
	Wum2_mat = spdiags(wum.^2,0,npix,npix);
	Wvp2_mat = spdiags(wvp.^2,0,npix,npix);
	Wvm2_mat = spdiags(wvm.^2,0,npix,npix);
	Eup2_mat = spdiags((Dup*z(imask)-p(imask)).^2,0,npix,npix);
	Eum2_mat = spdiags((Dum*z(imask)-p(imask)).^2,0,npix,npix);
	Evp2_mat = spdiags((Dvp*z(imask)-q(imask)).^2,0,npix,npix);
	Evm2_mat = spdiags((Dvm*z(imask)-q(imask)).^2,0,npix,npix);
	energie = norm(Lambda_two*(z(imask)-z0(imask)))+0.5*((norm(wup-1))/(4*epsilon)+(norm(wum-1))/(4*epsilon)+(norm(wvp-1))/(4*epsilon)+(norm(wvm-1))/(4*epsilon)+epsilon*norm(Dup*wup)+epsilon*norm(Dum*wum)+epsilon*norm(Dvp*wvp)+epsilon*norm(Dvm*wvp)+mu*norm(sqrt(Eup2_mat)*wup)+mu*norm(sqrt(Eum2_mat)*wum)+mu*norm(sqrt(Evp2_mat)*wvp)+mu*norm(sqrt(Evm2_mat)*wvm));
	if(nargout>5)
		tab_nrj(cpt) = energie;
		cpt = cpt+1;
	end
	if(nargout>6)
		lambda = -mean(z(imask)-gt(imask));
		zrmse = z+lambda;
		tab_rmse(cpt_rmse) = sqrt(mean((zrmse(imask)-gt(imask)).^2));
		cpt_rmse = cpt_rmse+1;
	end
	
	% Alternating optimisation loops
	for k = 1:maxit


		% w update
		Eup2_mat = spdiags((Dup*z(imask)-p(imask)).^2,0,npix,npix);
		Eum2_mat = spdiags((Dum*z(imask)-p(imask)).^2,0,npix,npix);
		Evp2_mat = spdiags((Dvp*z(imask)-q(imask)).^2,0,npix,npix);
		Evm2_mat = spdiags((Dvm*z(imask)-q(imask)).^2,0,npix,npix);
		Aup = 0.5*(mu*Eup2_mat+epsilon*DuptDup+(1/(4*epsilon))*speye(npix));
		Avp = 0.5*(mu*Evp2_mat+epsilon*DvptDvp+(1/(4*epsilon))*speye(npix));
		Aum = 0.5*(mu*Eum2_mat+epsilon*DumtDum+(1/(4*epsilon))*speye(npix));
		Avm = 0.5*(mu*Evm2_mat+epsilon*DvmtDvm+(1/(4*epsilon))*speye(npix));
		[wup,fl] = pcg(Aup,bw,1e-4,100,[],[],wup);
		[wum,fl] = pcg(Aum,bw,1e-4,100,[],[],wum);
		[wvp,fl] = pcg(Avp,bw,1e-4,100,[],[],wvp);
		[wvm,fl] = pcg(Avm,bw,1e-4,100,[],[],wvm);
		Wup2_mat = spdiags(wup.^2,0,npix,npix);
		Wum2_mat = spdiags(wum.^2,0,npix,npix);
		Wvp2_mat = spdiags(wvp.^2,0,npix,npix);
		Wvm2_mat = spdiags(wvm.^2,0,npix,npix);
				
		% z update

		A = 0.5*mu*(Dup'*Wup2_mat*Dup+Dum'*Wum2_mat*Dum+Dvp'*Wvp2_mat*Dvp+Dvm'*Wvm2_mat*Dvm)+Lambda_two; % Matrix of the system
		b =  0.5*mu*(Dup'*Wup2_mat+Dum'*Wum2_mat)*p(imask)+0.5*mu*(Dvp'*Wvp2_mat+Dvm'*Wvm2_mat)*q(imask)+Lambda_two*z0(imask);		
		[z(imask),fl] = pcg(A,b,1e-4,100,[],[],z(imask));
		


		% Check CV
		energie_old = energie;
		energie = norm(Lambda_two*(z(imask)-z0(imask)))+0.5*((norm(wup-1))/(4*epsilon)+(norm(wum-1))/(4*epsilon)+(norm(wvp-1))/(4*epsilon)+(norm(wvm-1))/(4*epsilon)+epsilon*norm(Dup*wup)+epsilon*norm(Dum*wum)+epsilon*norm(Dvp*wvp)+epsilon*norm(Dvm*wvp)+mu*norm(sqrt(Eup2_mat)*wup)+mu*norm(sqrt(Eum2_mat)*wum)+mu*norm(sqrt(Evp2_mat)*wvp)+mu*norm(sqrt(Evm2_mat)*wvm));
		if(nargout>5)
			tab_nrj(cpt) = energie;
			cpt = cpt+1;
		end
		if(nargout>6)
			lambda = -mean(z(imask)-gt(imask));
			zrmse = z+lambda;
			tab_rmse(cpt_rmse) = sqrt(mean((zrmse(imask)-gt(imask)).^2));
			cpt_rmse = cpt_rmse+1;
		end

			if(nargout>5)
			if(mod(k,20)==1)
			figure(3)
			plot(tab_nrj(1:cpt-1))
			%~ drawnow
			end
		end
		if(nargout>6)
		if(mod(k,20)==1)
			figure(4)
			plot(tab_rmse(1:cpt_rmse-1))
			%~ drawnow
		end
		end
		%~ if(mod(k,5)==1)
		%~ figure(1)
		%~ surfl(z,[-135 30]);
		%~ view(-35,20)		
		%~ shading flat;
		%~ colormap gray;
		%~ grid off
		%~ axis ij;
		%~ axis equal;
		%~ axis([1 size(p,1) 1 size(p,2) min(z(:)) max(z(:))]);
		%~ drawnow
		%~ end		
		
		relative_residual = abs(energie-energie_old)./abs(energie_old);
		disp(sprintf('it %d - EAT = %.4f - res : %.6f',k,energie,relative_residual))
		if(relative_residual < tol & k>5)
			break;
		end
	end
	if(k==maxit)
		disp('max number of iterations reached');
	end
	
	% Put NaNs outside the mask
	z(mask==0) = NaN;	
	wup(mask==0) = NaN;	
	wum(mask==0) = NaN;	
	wvp(mask==0) = NaN;	
	wvm(mask==0) = NaN;

	if(nargout>5)
		tab_nrj = tab_nrj(1:k+1);
	end
	if(nargout>6)
		tab_rmse = tab_rmse(1:k+1);
	end
end
