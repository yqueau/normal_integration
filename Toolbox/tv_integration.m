function [z,tab_nrj,tab_rmse] = tv_integration(p,q,mask,lambda,z0,alpha,maxit,tol,zinit,gt)

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
	if (~exist('alpha','var')|isempty(alpha)) alpha = 0.01; end;
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
	A = Lambda_two+0.125*alpha*((Dup'*Dup+Dvp'*Dvp)+(Dup'*Dup+Dvm'*Dvm)+(Dum'*Dum+Dvp'*Dvp)+(Dum'*Dum+Dvm'*Dvm));
	
	if(nargout>1)
		tab_nrj = zeros(maxit+1,1);
		cpt = 1;
	end
	if(nargout>2)
		tab_rmse = zeros(maxit+1,1);
		cpt_rmse = 1;
	end
	
	% Initialization
	z = zinit;

	rpp1 = zeros(npix,1);
	rpm1 = zeros(npix,1);
	rmp1 = zeros(npix,1);
	rmm1 = zeros(npix,1);
	rpp2 = zeros(npix,1);
	rpm2 = zeros(npix,1);
	rmp2 = zeros(npix,1);
	rmm2 = zeros(npix,1);
	bpp1 = zeros(npix,1);
	bpm1 = zeros(npix,1);
	bmp1 = zeros(npix,1);
	bmm1 = zeros(npix,1);
	bpp2 = zeros(npix,1);
	bpm2 = zeros(npix,1);
	bmp2 = zeros(npix,1);
	bmm2 = zeros(npix,1);

	% r update
	spp1 = Dup*z(imask)-p(imask)+bpp1; 
	spm1 = Dup*z(imask)-p(imask)+bpm1; 
	smp1 = Dum*z(imask)-p(imask)+bmp1; 
	smm1 = Dum*z(imask)-p(imask)+bmm1; 
	spp2 = Dvp*z(imask)-q(imask)+bpp2; 
	spm2 = Dvm*z(imask)-q(imask)+bpm2; 
	smp2 = Dvp*z(imask)-q(imask)+bmp2; 
	smm2 = Dvm*z(imask)-q(imask)+bmm2;
	rpp1 = max(sqrt(spp1.^2+spp2.^2)-4/alpha,0).*spp1./sqrt(spp1.^2+spp2.^2);
	rpm1 = max(sqrt(spm1.^2+spm2.^2)-4/alpha,0).*spm1./sqrt(spm1.^2+spm2.^2);
	rmp1 = max(sqrt(smp1.^2+smp2.^2)-4/alpha,0).*smp1./sqrt(smp1.^2+smp2.^2);
	rmm1 = max(sqrt(smm1.^2+smm2.^2)-4/alpha,0).*smm1./sqrt(smm1.^2+smm2.^2);
	rpp2 = max(sqrt(spp1.^2+spp2.^2)-4/alpha,0).*spp2./sqrt(spp1.^2+spp2.^2);
	rpm2 = max(sqrt(spm1.^2+spm2.^2)-4/alpha,0).*spm2./sqrt(spm1.^2+spm2.^2);
	rmp2 = max(sqrt(smp1.^2+smp2.^2)-4/alpha,0).*smp2./sqrt(smp1.^2+smp2.^2);
	rmm2 = max(sqrt(smm1.^2+smm2.^2)-4/alpha,0).*smm2./sqrt(smm1.^2+smm2.^2);

	% b update
	bpp1 = bpp1+Dup*z(imask)-p(imask)-rpp1;
	bpm1 = bpm1+Dup*z(imask)-p(imask)-rpm1;
	bmp1 = bmp1+Dum*z(imask)-p(imask)-rmp1;
	bmm1 = bmm1+Dum*z(imask)-p(imask)-rmm1;
	bpp2 = bpp2+Dvp*z(imask)-q(imask)-rpp2;
	bpm2 = bpm2+Dvm*z(imask)-q(imask)-rpm2;
	bmp2 = bmp2+Dvp*z(imask)-q(imask)-rmp2;
	bmm2 = bmm2+Dvm*z(imask)-q(imask)-rmm2;	
	
	energie = norm(Lambda_two*(z(imask)-z0(imask)))+sum(sqrt((Dup*z(imask)-p(imask)).^2)+(Dvp*z(imask)-q(imask)).^2)+sum(sqrt((Dup*z(imask)-p(imask)).^2)+(Dvm*z(imask)-q(imask)).^2)+sum(sqrt((Dum*z(imask)-p(imask)).^2)+(Dvp*z(imask)-q(imask)).^2)+sum(sqrt((Dum*z(imask)-p(imask)).^2)+(Dvm*z(imask)-q(imask)).^2);
	
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
		ppp = p(imask)+rpp1-bpp1;
		ppm = p(imask)+rpm1-bpm1;
		pmp = p(imask)+rmp1-bmp1;
		pmm = p(imask)+rmm1-bmm1;
		qpp = q(imask)+rpp2-bpp2;
		qpm = q(imask)+rpm2-bpm2;
		qmp = q(imask)+rmp2-bmp2;
		qmm = q(imask)+rmm2-bmm2;
		b = Lambda_two*z0(imask)+0.125*alpha*((Dup'*ppp+Dvp'*qpp)+(Dup'*ppm+Dvm'*qpm)+(Dum'*pmp+Dvp'*qmp)+(Dum'*pmm+Dvm'*qmm));
		[z(imask),fl] = pcg(A,b,1e-9,10,[],[],z(imask));

		% r update
		spp1 = Dup*z(imask)-p(imask)+bpp1; 
		spm1 = Dup*z(imask)-p(imask)+bpm1; 
		smp1 = Dum*z(imask)-p(imask)+bmp1; 
		smm1 = Dum*z(imask)-p(imask)+bmm1; 
		spp2 = Dvp*z(imask)-q(imask)+bpp2; 
		spm2 = Dvm*z(imask)-q(imask)+bpm2; 
		smp2 = Dvp*z(imask)-q(imask)+bmp2; 
		smm2 = Dvm*z(imask)-q(imask)+bmm2;
		rpp1 = max(sqrt(spp1.^2+spp2.^2)-4/alpha,0).*spp1./sqrt(spp1.^2+spp2.^2);
		rpm1 = max(sqrt(spm1.^2+spm2.^2)-4/alpha,0).*spm1./sqrt(spm1.^2+spm2.^2);
		rmp1 = max(sqrt(smp1.^2+smp2.^2)-4/alpha,0).*smp1./sqrt(smp1.^2+smp2.^2);
		rmm1 = max(sqrt(smm1.^2+smm2.^2)-4/alpha,0).*smm1./sqrt(smm1.^2+smm2.^2);
		rpp2 = max(sqrt(spp1.^2+spp2.^2)-4/alpha,0).*spp2./sqrt(spp1.^2+spp2.^2);
		rpm2 = max(sqrt(spm1.^2+spm2.^2)-4/alpha,0).*spm2./sqrt(spm1.^2+spm2.^2);
		rmp2 = max(sqrt(smp1.^2+smp2.^2)-4/alpha,0).*smp2./sqrt(smp1.^2+smp2.^2);
		rmm2 = max(sqrt(smm1.^2+smm2.^2)-4/alpha,0).*smm2./sqrt(smm1.^2+smm2.^2);

		% b update
		bpp1 = bpp1+Dup*z(imask)-p(imask)-rpp1;
		bpm1 = bpm1+Dup*z(imask)-p(imask)-rpm1;
		bmp1 = bmp1+Dum*z(imask)-p(imask)-rmp1;
		bmm1 = bmm1+Dum*z(imask)-p(imask)-rmm1;
		bpp2 = bpp2+Dvp*z(imask)-q(imask)-rpp2;
		bpm2 = bpm2+Dvm*z(imask)-q(imask)-rpm2;
		bmp2 = bmp2+Dvp*z(imask)-q(imask)-rmp2;
		bmm2 = bmm2+Dvm*z(imask)-q(imask)-rmm2;

		
		% Check CV
		energie_old = energie;
		energie = norm(Lambda_two*(z(imask)-z0(imask)))+0.25*(sum(sqrt((Dup*z(imask)-p(imask)).^2)+(Dvp*z(imask)-q(imask)).^2)+sum(sqrt((Dup*z(imask)-p(imask)).^2)+(Dvm*z(imask)-q(imask)).^2)+sum(sqrt((Dum*z(imask)-p(imask)).^2)+(Dvp*z(imask)-q(imask)).^2)+sum(sqrt((Dum*z(imask)-p(imask)).^2)+(Dvm*z(imask)-q(imask)).^2));
	
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
%~ 
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
		%~ wupdisp(imask) = bpp1;
		%~ imagesc(wupdisp);
		%~ colormap gray
		%~ colorbar
%~ 
		%~ if(nargout>1)
			%~ figure(3)
			%~ plot(tab_nrj(1:cpt-1))
			%~ drawnow
		%~ end
		%~ if(nargout>2)
			%~ figure(4)
			%~ plot(tab_rmse(1:cpt_rmse-1))
			%~ drawnow
		%~ end		
		
		relative_residual = abs(energie-energie_old)./abs(energie_old);
		disp(sprintf('it %d - EAT = %.4f - res : %.6f',k,energie,relative_residual));
		if(relative_residual < tol & k>20)
			break;
		end
	end

	if(k==maxit)
		disp('Maximum number of iterations reached');
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
