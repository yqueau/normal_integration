function [z,tab_nrj,tab_rmse] = phi1_integration(p,q,mask,lambda,z0,beta,gamma,maxit,tol,zinit,gt)

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
	if (~exist('beta','var')|isempty(beta)) beta = 0.99; end;
	if (~exist('gamma','var')|isempty(gamma)) gamma = 0.01; end;
	if (~exist('maxit','var')|isempty(maxit)) maxit = 100; end;
	if (~exist('tol','var')|isempty(tol)) tol = 1e-4; end;
	if (~exist('zinit','var')|isempty(zinit)) zinit = z0; end;
	if (~exist('gt','var')|isempty(gt)) gt = z0; end;
	

	max_backtracking = 100;
	Linit = 0.5;
	L = Linit;
	c = 1e-2;
	eta = 1.2;

	% If lambda is a scalar, make it a matrix
	if(size(lambda,1)==1)
		lambda = lambda*mask;
	end

	% Make finite differences operators	
	[Dup,Dum,Dvp,Dvm,imask,imaskup,imaskum,imaskvp,imaskvm] = make_gradient(mask);
	npix = length(imask);

	% Some stuff used later
	II = transpose(1:2*npix);
	JJ = repmat(1:npix,[2 1]);
	JJ = JJ(:); % 1 1 2 2 3 3 4 4....
	Dt = sparse([],[],[],npix,2*npix,6*npix);
	Dt(:,1:2:end-1) = Dup';
	Dt(:,2:2:end) = Dvp';
	Dt2 = sparse([],[],[],npix,2*npix,6*npix);
	Dt2(:,1:2:end-1) = Dup';
	Dt2(:,2:2:end) = Dvm';
	Dt3 = sparse([],[],[],npix,2*npix,6*npix);
	Dt3(:,1:2:end-1) = Dum';
	Dt3(:,2:2:end) = Dvp';
	Dt4 = sparse([],[],[],npix,2*npix,6*npix);
	Dt4(:,1:2:end-1) = Dum';
	Dt4(:,2:2:end) = Dvm';
	Dz_minus_g = zeros(2,npix);
	Dz_minus_g2 = zeros(2,npix);
	Dz_minus_g3 = zeros(2,npix);
	Dz_minus_g4 = zeros(2,npix);

	Lambda_two = spdiags(lambda(imask),0,npix,npix); % Regularization

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
	zprevious = z;	


	% Compute gradient of f at current estimate
	% Residual
	Dz_minus_g(1,:) = Dup*z(imask)-p(imask);
	Dz_minus_g(2,:) = Dvp*z(imask)-q(imask);
	Dz_minus_g2(1,:) = Dup*z(imask)-p(imask);
	Dz_minus_g2(2,:) = Dvm*z(imask)-q(imask);
	Dz_minus_g3(1,:) = Dum*z(imask)-p(imask);
	Dz_minus_g3(2,:) = Dvp*z(imask)-q(imask);
	Dz_minus_g4(1,:) = Dum*z(imask)-p(imask);
	Dz_minus_g4(2,:) = Dvm*z(imask)-q(imask);
	% Normalized residual
	norme_Dz_minus_g = gamma^2+(Dz_minus_g(1,:).^2+Dz_minus_g(2,:).^2);
	norme_Dz_minus_g2 = gamma^2+(Dz_minus_g2(1,:).^2+Dz_minus_g2(2,:).^2);
	norme_Dz_minus_g3 = gamma^2+(Dz_minus_g3(1,:).^2+Dz_minus_g3(2,:).^2);
	norme_Dz_minus_g4 = gamma^2+(Dz_minus_g4(1,:).^2+Dz_minus_g4(2,:).^2);

	energie = norm(Lambda_two*(z(imask)-z0(imask)))+0.25*(sum(log(norme_Dz_minus_g))+sum(log(norme_Dz_minus_g2))+sum(log(norme_Dz_minus_g3))+sum(log(norme_Dz_minus_g4)));
	
	if(nargout>1)
		tab_nrj(cpt) = energie;
		cpt = cpt+1;
	end
	if(nargout>2)
		lambda_rmse = -mean(z(imask)-gt(imask));
		zrmse = z+lambda_rmse;
		tab_rmse(cpt_rmse) = sqrt(mean((zrmse(imask)-gt(imask)).^2));
		cpt_rmse = cpt_rmse+1;
	end

	% Current energy
	f_curr = 0.25*(sum(log(norme_Dz_minus_g))+sum(log(norme_Dz_minus_g2))+sum(log(norme_Dz_minus_g3))+sum(log(norme_Dz_minus_g4)));

	% Alternating optimisation loops
	for k = 1:maxit

		Lcurr = L;

		% Current gradient
		Dz_minus_g = bsxfun(@rdivide,Dz_minus_g,norme_Dz_minus_g);
		Dz_minus_g2 = bsxfun(@rdivide,Dz_minus_g2,norme_Dz_minus_g2);
		Dz_minus_g3 = bsxfun(@rdivide,Dz_minus_g3,norme_Dz_minus_g3);
		Dz_minus_g4 = bsxfun(@rdivide,Dz_minus_g4,norme_Dz_minus_g4);
		% Make it a matrix
		Dz_minus_g_mat = sparse(II,JJ,Dz_minus_g(:));
		Dz_minus_g_mat2 = sparse(II,JJ,Dz_minus_g2(:));
		Dz_minus_g_mat3 = sparse(II,JJ,Dz_minus_g3(:));
		Dz_minus_g_mat4 = sparse(II,JJ,Dz_minus_g4(:));
		% Get all terms inside the sum
		grad_f_curr = sum(Dt*Dz_minus_g_mat+Dt2*Dz_minus_g_mat2+Dt3*Dz_minus_g_mat3+Dt4*Dz_minus_g_mat4,2);

		% Lazy backtracking to set stepsize
		lc = 0; % lc = 1 if Lipschitz constant L is big enough
		while(lc < max_backtracking)
			alpha = 1.99*(1-beta)/L; % Descent stepsize

			% Forward update, given the gradient
			zbar = z(imask) - alpha*grad_f_curr+beta*(z(imask)-zprevious(imask));
			% Backward update (prox. update)
			znext = (zbar + 2*alpha*(lambda(imask)).*z0(imask))./(1+2*alpha*(lambda(imask)));

			z_dist = znext-z(imask); % Evaluate the difference between current and next estimate

			% Next energy
			Dz_minus_g(1,:) = Dup*znext(imask)-p(imask);
			Dz_minus_g(2,:) = Dvp*znext(imask)-q(imask);
			Dz_minus_g2(1,:) = Dup*znext(imask)-p(imask);
			Dz_minus_g2(2,:) = Dvm*znext(imask)-q(imask);
			Dz_minus_g3(1,:) = Dum*znext(imask)-p(imask);
			Dz_minus_g3(2,:) = Dvp*znext(imask)-q(imask);
			Dz_minus_g4(1,:) = Dum*znext(imask)-p(imask);
			Dz_minus_g4(2,:) = Dvm*znext(imask)-q(imask);
			norme_Dz_minus_g = gamma^2+(Dz_minus_g(1,:).^2+Dz_minus_g(2,:).^2);
			norme_Dz_minus_g2 = gamma^2+(Dz_minus_g2(1,:).^2+Dz_minus_g2(2,:).^2);
			norme_Dz_minus_g3 = gamma^2+(Dz_minus_g3(1,:).^2+Dz_minus_g3(2,:).^2);
			norme_Dz_minus_g4 = gamma^2+(Dz_minus_g4(1,:).^2+Dz_minus_g4(2,:).^2);

			f_next = 0.25*(sum(log(norme_Dz_minus_g))+sum(log(norme_Dz_minus_g2))+sum(log(norme_Dz_minus_g3))+sum(log(norme_Dz_minus_g4)));				
			
			% Lipschitz test
			if(f_next <= f_curr+grad_f_curr'*z_dist+0.5*L*z_dist'*z_dist)
				L = L/1.05;
				break; % if Lipschitz => stepsize is small enough
			else
				lc = lc+1; % if not Lipschitz => try smaller stepsize
				L = eta*L;
			end
		end

		% Update auxiliary variables
		zprevious = z;
		z(imask) = znext;
		f_curr = f_next;
		
		% Current residual
		Dz_minus_g(1,:) = Dup*z(imask)-p(imask);
		Dz_minus_g(2,:) = Dvp*z(imask)-q(imask);
		Dz_minus_g2(1,:) = Dup*z(imask)-p(imask);
		Dz_minus_g2(2,:) = Dvm*z(imask)-q(imask);
		Dz_minus_g3(1,:) = Dum*z(imask)-p(imask);
		Dz_minus_g3(2,:) = Dvp*z(imask)-q(imask);
		Dz_minus_g4(1,:) = Dum*z(imask)-p(imask);
		Dz_minus_g4(2,:) = Dvm*z(imask)-q(imask);
		% Normalized residual
		norme_Dz_minus_g = gamma^2+(Dz_minus_g(1,:).^2+Dz_minus_g(2,:).^2);
		norme_Dz_minus_g2 = gamma^2+(Dz_minus_g2(1,:).^2+Dz_minus_g2(2,:).^2);
		norme_Dz_minus_g3 = gamma^2+(Dz_minus_g3(1,:).^2+Dz_minus_g3(2,:).^2);
		norme_Dz_minus_g4 = gamma^2+(Dz_minus_g4(1,:).^2+Dz_minus_g4(2,:).^2);
			
%~ 
		%~ % After a few iterations, decrease the Lipschitz constant for speedup
		%~ if(mod(k,50) == 1)
%~ 
			%~ % Current gradient
			%~ Dz_minus_g = bsxfun(@rdivide,Dz_minus_g,norme_Dz_minus_g);
			%~ Dz_minus_g2 = bsxfun(@rdivide,Dz_minus_g2,norme_Dz_minus_g2);
			%~ Dz_minus_g3 = bsxfun(@rdivide,Dz_minus_g3,norme_Dz_minus_g3);
			%~ Dz_minus_g4 = bsxfun(@rdivide,Dz_minus_g4,norme_Dz_minus_g4);
			%~ % Make it a matrix
			%~ Dz_minus_g_mat = sparse(II,JJ,Dz_minus_g(:));
			%~ Dz_minus_g_mat2 = sparse(II,JJ,Dz_minus_g2(:));
			%~ Dz_minus_g_mat3 = sparse(II,JJ,Dz_minus_g3(:));
			%~ Dz_minus_g_mat4 = sparse(II,JJ,Dz_minus_g4(:));
			%~ % Get all terms inside the sum
			%~ grad_f_curr = sum(Dt*Dz_minus_g_mat+Dt2*Dz_minus_g_mat2+Dt3*Dz_minus_g_mat3+Dt4*Dz_minus_g_mat4,2);
			%~ 
			%~ % Lazy backtracking to set stepsize
			%~ lc = 0; % lc = 1 if Lipschitz constant L is big enough
			%~ while(lc < max_backtracking)
				%~ alpha = 2*(1-beta)/(c+L); % Descent stepsize
%~ 
				%~ % Forward update, given the gradient
				%~ zbar = z(imask) - alpha*grad_f_curr+beta*(z(imask)-zprevious(imask));
				%~ % Backward update (prox. update)
				%~ znext = (zbar + 2*alpha*(lambda(imask)).*z0(imask))./(1+2*alpha*(lambda(imask)));
%~ 
				%~ z_dist = znext-z(imask); % Evaluate the difference between current and next estimate
%~ 
				%~ % Next energy
				%~ Dz_minus_g(1,:) = Dup*znext(imask)-p(imask);
				%~ Dz_minus_g(2,:) = Dvp*znext(imask)-q(imask);
				%~ Dz_minus_g2(1,:) = Dup*znext(imask)-p(imask);
				%~ Dz_minus_g2(2,:) = Dvm*znext(imask)-q(imask);
				%~ Dz_minus_g3(1,:) = Dum*znext(imask)-p(imask);
				%~ Dz_minus_g3(2,:) = Dvp*znext(imask)-q(imask);
				%~ Dz_minus_g4(1,:) = Dum*znext(imask)-p(imask);
				%~ Dz_minus_g4(2,:) = Dvm*znext(imask)-q(imask);
				%~ norme_Dz_minus_g = gamma^2+(Dz_minus_g(1,:).^2+Dz_minus_g(2,:).^2);
				%~ norme_Dz_minus_g2 = gamma^2+(Dz_minus_g2(1,:).^2+Dz_minus_g2(2,:).^2);
				%~ norme_Dz_minus_g3 = gamma^2+(Dz_minus_g3(1,:).^2+Dz_minus_g3(2,:).^2);
				%~ norme_Dz_minus_g4 = gamma^2+(Dz_minus_g4(1,:).^2+Dz_minus_g4(2,:).^2);
%~ 
				%~ f_next = 0.5*(sum(log(norme_Dz_minus_g))+sum(log(norme_Dz_minus_g2))+sum(log(norme_Dz_minus_g3))+sum(log(norme_Dz_minus_g4)));				
				%~ 
				%~ % Lipschitz test
				%~ if(f_next <= f_curr+grad_f_curr'*z_dist+0.5*L*z_dist'*z_dist)
					%~ L = L/eta;
					%~ lc = lc+1; % if not Lipschitz => try smaller stepsize					
				%~ else					
					%~ lc = max_backtracking; % if Lipschitz => stepsize is small enough					
				%~ end
			%~ end
		%~ end


		% Check CV
		energie_old = energie;
		energie = norm(Lambda_two*(z(imask)-z0(imask)))+0.25*(sum(log(norme_Dz_minus_g))+sum(log(norme_Dz_minus_g2))+sum(log(norme_Dz_minus_g3))+sum(log(norme_Dz_minus_g4)));
	
		if(nargout>1)
			tab_nrj(cpt) = energie;
			cpt = cpt+1;
		end
		if(nargout>2)
			lambda_rmse = -mean(z(imask)-gt(imask));
			zrmse = z+lambda_rmse;
			tab_rmse(cpt_rmse) = sqrt(mean((zrmse(imask)-gt(imask)).^2));
			cpt_rmse = cpt_rmse+1;
		end


		if(nargout>1)
			if(mod(k,100)==1)
			figure(3)
			plot(tab_nrj(1:cpt-1))
			%~ drawnow
			end
		end
		if(nargout>2)
		if(mod(k,100)==1)
			figure(4)
			plot(tab_rmse(1:cpt_rmse-1))
			%~ drawnow
		end
		end
		%~ if(mod(k,50)==1)
		%~ figure(473)
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
		disp(sprintf('it %d - EAT = %.4f - res : %.6f',k,energie,relative_residual));
		if(relative_residual < tol & (energie<energie_old) & k > 500)
			break;
		end
	end
	if(k == maxit)
		disp('Max number of iterations reached');
	end

	%~ close(473)
	
	% Put NaNs outside the mask
	z(mask==0) = NaN;	

	if(nargout>1)
		tab_nrj = tab_nrj(1:k+1);
	end
	if(nargout>2)
		tab_rmse = tab_rmse(1:k+1);
	end
end
