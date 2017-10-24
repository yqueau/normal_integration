function z = smooth_integration(p,q,mask,lambda,z0,solver,precond)

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
	if (~exist('solver','var')|isempty(solver)) solver = 'pcg'; end;
	if (~exist('precond','var')|isempty(precond)) precond = 'ichol'; end;

	% If lambda is a scalar, make it a matrix
	if(size(lambda,1)==1)
		lambda = lambda*mask;
	end

	% Make finite differences operators	
	[Dup,Dum,Dvp,Dvm,imask] = make_gradient(mask);

	% Matrix of the system
	L = 0.5*(Dup'*Dup+Dum'*Dum+Dvp'*Dvp+Dvm'*Dvm);
	Lambda_two = spdiags(lambda(imask),0,length(imask),length(imask)); 
	A = L + Lambda_two;

	% Second membre
	Du = 0.5*(Dup'+Dum');
	Dv = 0.5*(Dvp'+Dvm');
	b = Du*p(imask)+Dv*q(imask)+Lambda_two*z0(imask);

	% Preconditioning
	if(strcmp(precond,'none'))
		precondL = [];
		precondR = [];
	elseif(strcmp(precond,'CMG'))
		precondL = cmg_sdd(A);
		precondR = [];
	elseif(strcmp(precond,'ichol')) % Modified incomplete cholesky advised in [Bahr et al., CVM 2017]
		precondL = ichol(A,struct('type','ict','droptol',1e-03,'michol','on'));	
		precondR = precondL';	
	end

	% Resolution
	z = z0;
	if(strcmp(solver,'direct')) % Calls cholesky
		z(imask) = A\b;
	elseif(strcmp(solver,'pcg')) % Calls CG
		z(imask) = pcg(A,b,1e-4,1000,precondL,precondR,z(imask));
	end

	% Put NaNs outside the mask
	z(mask==0) = NaN;	
end
