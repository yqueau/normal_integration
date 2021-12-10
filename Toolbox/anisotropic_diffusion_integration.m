function [z,tab_nrj,tab_rmse] = anisotropic_diffusion_integration(p,q,mask,lambda,z0,mu,nu,maxit,tol,zinit,gt)

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
	if (~exist('mu','var')|isempty(mu)) mu = 0.01; end;
	if (~exist('nu','var')|isempty(nu)) nu = 0.01; end;
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
	Lambda_two = spdiags(lambda(imask),0,npix,npix); % Regularization
	if(nargout>1)
		tab_nrj = zeros(maxit+1,1);
		cpt = 1;
	end
	if(nargout>2)
		tab_rmse = zeros(maxit+1,1);
		cpt_rmse = 1;
	end
	nu2 = nu*nu;

	% Initialization
	z = zinit;
	App_mat = spdiags(1./(sqrt(1+p(imask).^2./nu2).*sqrt(1+(Dup*z(imask)./mu).^2+(Dvp*z(imask)./mu).^2)),0,npix,npix);
	Apm_mat = spdiags(1./(sqrt(1+p(imask).^2./nu2).*sqrt(1+(Dup*z(imask)./mu).^2+(Dvm*z(imask)./mu).^2)),0,npix,npix);
	Amp_mat = spdiags(1./(sqrt(1+p(imask).^2./nu2).*sqrt(1+(Dum*z(imask)./mu).^2+(Dvp*z(imask)./mu).^2)),0,npix,npix);
	Amm_mat = spdiags(1./(sqrt(1+p(imask).^2./nu2).*sqrt(1+(Dum*z(imask)./mu).^2+(Dvm*z(imask)./mu).^2)),0,npix,npix);
	Bpp_mat = spdiags(1./(sqrt(1+q(imask).^2./nu2).*sqrt(1+(Dup*z(imask)./mu).^2+(Dvp*z(imask)./mu).^2)),0,npix,npix);
	Bpm_mat = spdiags(1./(sqrt(1+q(imask).^2./nu2).*sqrt(1+(Dup*z(imask)./mu).^2+(Dvm*z(imask)./mu).^2)),0,npix,npix);
	Bmp_mat = spdiags(1./(sqrt(1+q(imask).^2./nu2).*sqrt(1+(Dum*z(imask)./mu).^2+(Dvp*z(imask)./mu).^2)),0,npix,npix);
	Bmm_mat = spdiags(1./(sqrt(1+q(imask).^2./nu2).*sqrt(1+(Dum*z(imask)./mu).^2+(Dvm*z(imask)./mu).^2)),0,npix,npix);
	energie = norm(Lambda_two*(z(imask)-z0(imask)))+norm(App_mat*(Dup*z(imask)-p(imask)))+norm(Bpp_mat*(Dvp*z(imask)-q(imask)))+norm(Apm_mat*(Dup*z(imask)-p(imask)))+norm(Bpm_mat*(Dvm*z(imask)-q(imask)))+norm(Amp_mat*(Dum*z(imask)-p(imask)))+norm(Bmp_mat*(Dvp*z(imask)-q(imask)))+norm(Amm_mat*(Dum*z(imask)-p(imask)))+norm(Bmm_mat*(Dvm*z(imask)-q(imask)));
	if(nargout>1)
		tab_nrj(cpt) = energie;
		cpt = cpt+1;
	end
	if(nargout>2)
		lambda = -mean(z(imask)-gt(imask));
		zrmse = z+lambda;
		tab_rmse(cpt_rmse) = sqrt(mean((zrmse(imask)-gt(imask)).^2));
		cpt_rmse = cpt_rmse+1;
	end
	
	% Alternating optimisation loops
	for k = 1:maxit
		
		% z update
		A = Lambda_two+(Dup'*(App_mat.^2)*Dup)+(Dvp'*(Bpp_mat.^2)*Dvp)+(Dup'*(Apm_mat.^2)*Dup)+(Dvm'*(Bpm_mat.^2)*Dvm)+(Dum'*(Amp_mat.^2)*Dum)+(Dvp'*(Bmp_mat.^2)*Dvp)+(Dum'*(Amm_mat.^2)*Dum)+(Dvm'*(Bmm_mat.^2)*Dvm);
		b = Lambda_two*z0(imask)+(Dup'*(App_mat.^2)*p(imask))+(Dvp'*(Bpp_mat.^2)*q(imask))+(Dup'*(Apm_mat.^2)*p(imask))+(Dvm'*(Bpm_mat.^2)*q(imask))+(Dum'*(Amp_mat.^2)*p(imask))+(Dvp'*(Bmp_mat.^2)*q(imask))+(Dum'*(Amm_mat.^2)*p(imask))+(Dvm'*(Bmm_mat.^2)*q(imask));
		%~ precond = cmg_sdd(A);
		%~ [z(imask)] = pcg(A,b,1e-4,100,precond,[],z(imask));
		[z(imask)] = A\b;
		
		% w update
		App_mat = spdiags(1./(sqrt(1+p(imask).^2./nu2).*sqrt(1./mu+(Dup*z(imask)./mu).^2+(Dvp*z(imask)./mu).^2)),0,npix,npix);
		Apm_mat = spdiags(1./(sqrt(1+p(imask).^2./nu2).*sqrt(1./mu+(Dup*z(imask)./mu).^2+(Dvm*z(imask)./mu).^2)),0,npix,npix);
		Amp_mat = spdiags(1./(sqrt(1+p(imask).^2./nu2).*sqrt(1./mu+(Dum*z(imask)./mu).^2+(Dvp*z(imask)./mu).^2)),0,npix,npix);
		Amm_mat = spdiags(1./(sqrt(1+p(imask).^2./nu2).*sqrt(1./mu+(Dum*z(imask)./mu).^2+(Dvm*z(imask)./mu).^2)),0,npix,npix);
		Bpp_mat = spdiags(1./(sqrt(1+q(imask).^2./nu2).*sqrt(1./mu+(Dup*z(imask)./mu).^2+(Dvp*z(imask)./mu).^2)),0,npix,npix);
		Bpm_mat = spdiags(1./(sqrt(1+q(imask).^2./nu2).*sqrt(1./mu+(Dup*z(imask)./mu).^2+(Dvm*z(imask)./mu).^2)),0,npix,npix);
		Bmp_mat = spdiags(1./(sqrt(1+q(imask).^2./nu2).*sqrt(1./mu+(Dum*z(imask)./mu).^2+(Dvp*z(imask)./mu).^2)),0,npix,npix);
		Bmm_mat = spdiags(1./(sqrt(1+q(imask).^2./nu2).*sqrt(1./mu+(Dum*z(imask)./mu).^2+(Dvm*z(imask)./mu).^2)),0,npix,npix);
		
		% Check CV
		energie_old = energie;
		energie = norm(Lambda_two*(z(imask)-z0(imask)))+norm(App_mat*(Dup*z(imask)-p(imask)))+norm(Bpp_mat*(Dvp*z(imask)-q(imask)))+norm(Apm_mat*(Dup*z(imask)-p(imask)))+norm(Bpm_mat*(Dvm*z(imask)-q(imask)))+norm(Amp_mat*(Dum*z(imask)-p(imask)))+norm(Bmp_mat*(Dvp*z(imask)-q(imask)))+norm(Amm_mat*(Dum*z(imask)-p(imask)))+norm(Bmm_mat*(Dvm*z(imask)-q(imask)));
		if(nargout>1)
			tab_nrj(cpt) = energie;
			cpt = cpt+1;
		end
		if(nargout>2)
			lambda = -mean(z(imask)-gt(imask));
			zrmse = z+lambda;
			tab_rmse(cpt_rmse) = sqrt(mean((zrmse(imask)-gt(imask)).^2));
			cpt_rmse = cpt_rmse+1;
		end

		%~ figure(1)
		%~ surfl(z,[-135 30]);
		%~ view(-35,20)		
		%~ axis ij;
		%~ axis equal;
		%~ axis([1 320 1 320 -20 70]);
		%~ shading flat;
		%~ colormap gray;
		%~ 
		%~ grid off
		%~ drawnow
%~ 
		%~ figure(2)
		%~ wupdisp = NaN*ones(nrows,ncols);
		%~ wupdisp(imask) = spdiags(App_mat,0);
		%~ imagesc(wupdisp);
		%~ colormap gray
		%~ colorbar

		if(nargout>1)
			figure(3)
			plot(tab_nrj(1:cpt-1))
			drawnow
		end
		if(nargout>2)
			figure(4)
			plot(tab_rmse(1:cpt_rmse-1))
			drawnow
		end		
		
		relative_residual = abs(energie-energie_old)./abs(energie_old);
		disp(sprintf('it %d - EAT = %.4f - res : %.6f',k,energie,relative_residual));
		if(relative_residual < tol & k > 5)	
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

	if(nargout>1)
		tab_nrj = tab_nrj(1:k+1);
	end
	if(nargout>2)
		tab_rmse = tab_rmse(1:k+1);
	end
end
