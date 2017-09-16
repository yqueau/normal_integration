function [u,ma,tab_rmse] = horn_brooks(p,q,omega,it_max,tol,trace,u0,ground_truth)
%hb integrates the gradient field [p,q] by minimizing the 
%functional F=\iint_{\Omega} \| \nabla U(x,y) - [p,q](x,y) \|^2 dx dy
%using the improved Horn and Brooks scheme
%	
%	U = hb(P,Q) uses default values
%	[U,ma,tab_rmse] = hb(P,Q,[],[],[],[],GT) also provides the computed masks and the evolution of RMSE between U and GT (ground truth)
%	U = hb(P,Q,OMEGA) uses the integration domain OMEGA (default : ones(size(p)))
%	U = hb(P,Q,OMEGA) uses the value gamma=ALPHA (default : 1)
%	U = hb(P,Q,OMEGA,IT_MAX) performs IT_MAX iterations (default : 100*size(p,1))
%	U = hb(P,Q,OMEGA,IT_MAX,TRACE) : if trace=1, it shows the recovered surface every 100 iterations (default : 0)
%	U = hb(P,Q,OMEGA,IT_MAX,TRACE,U0) : starts with surface U=U0 (default : zeros(size(p)))
%	U = hb(P,Q,OMEGA,IT_MAX,TRACE,U0,GT) : if a ground-truth is provided, the RMSE can be computed and displayed if TRACE=1			
	if (~exist('ground_truth','var')|isempty(ground_truth))
		save_rmse=0;
	else
		save_rmse=1;
	end
	
	if (~exist('tol','var')|isempty(tol)) tol=1e-3; end;
	if (~exist('u0','var')|isempty(u0)) u0=zeros(size(p)); end;
	if (~exist('trace','var')|isempty(trace)) trace=0; end;
	if (~exist('it_max','var')|isempty(it_max)) it_max=100*size(p,1); end;
	if (~exist('omega','var')|isempty(omega)) omega=ones(size(p)); end;
	
		
	

	[nrows,ncols]=size(omega);
	imask = find(omega>0);
	npix=length(imask);
	u=u0;
	u(omega==0)=NaN;
	tab_rmse=[];
	
	% Calcul des masques : 
	% ma1 : voisins de dessous et de droite dans omega
	% ma2 : ma0 \ ma1
	% ma3 : points de ma1 avec voisins de dessus et de gauche dans ma1
	% ma4 :                                      dans ma1, mais pas celui de gauche
	% ma5 :                               gauche dans ma1, mais pas celui de dessus
	% ma6 :                    ni voisin de gauche, ni voisin de dessus dans ma1
	% ma7 :           ma2 avec voisins de dessus et de gauche dans ma1
	% ma8 :                                      dans ma1, mais pas celui de gauche
	% ma9 :                               gauche dans ma1, mais pas celui de dessus
	% ma10:                     ni voisin de gauche, ni voisin de dessus dans ma1
	ma=zeros(nrows,ncols,9);
	ma(1:end-1,1:end-1,1)=omega(1:end-1,1:end-1).*omega(2:end,1:end-1).*omega(1:end-1,2:end);
	ma(:,:,2)=omega.*(~ma(:,:,1));
	ma(2:end,2:end,3)=ma(2:end,2:end,1).*ma(1:end-1,2:end,1).*ma(2:end,1:end-1,1);	
	ma(2:end,2:end,4)=ma(2:end,2:end,1).*ma(1:end-1,2:end,1).*(~ma(2:end,1:end-1,1));	
	ma(2:end,1,4)=ma(2:end,1,1).*ma(1:end-1,1,1);
	ma(2:end,2:end,5)=ma(2:end,2:end,1).*(~ma(1:end-1,2:end,1)).*(ma(2:end,1:end-1,1));	
	ma(1,2:end,5)=ma(1,2:end,1).*ma(1,1:end-1,1);
	ma(2:end,2:end,6)=ma(2:end,2:end,1).*(~ma(1:end-1,2:end,1)).*(~ma(2:end,1:end-1,1));	
	ma(2:end,1,6)=ma(2:end,1,1).*(~ma(1:end-1,1,1)); % Added	
	ma(1,2:end,6)=ma(1,2:end,1).*(~ma(1,1:end-1,1)); % Added	
	ma(1,1,6)=ma(1,1,1);
	ma(2:end,2:end,7)=ma(2:end,2:end,2).*(ma(1:end-1,2:end,1)).*(ma(2:end,1:end-1,1));	
	ma(2:end,2:end,8)=ma(2:end,2:end,2).*(ma(1:end-1,2:end,1)).*(~ma(2:end,1:end-1,1));	
	ma(2:end,1,8)=ma(2:end,1,2).*ma(1:end-1,1,1);
	ma(2:end,2:end,9)=ma(2:end,2:end,2).*(~ma(1:end-1,2:end,1)).*(ma(2:end,1:end-1,1));		
	ma(1,2:end,9)=ma(1,2:end,2).*ma(1,1:end-1,1);
	
	ind3=find(ma(:,:,3)>0);
	ind4=find(ma(:,:,4)>0);
	ind5=find(ma(:,:,5)>0);
	ind6=find(ma(:,:,6)>0);
	ind7=find(ma(:,:,7)>0);
	ind8=find(ma(:,:,8)>0);
	ind9=find(ma(:,:,9)>0);
	
	p_haut=[zeros(1,ncols);p(1:end-1,:)];
	p_bas =[p(2:end,:);zeros(1,ncols)];
	q_gauche=[zeros(nrows,1),q(:,1:end-1)];
	q_droite=[q(:,2:end),zeros(nrows,1)];	
	pq3=0.125*(p_haut(ind3)-p_bas(ind3)+q_gauche(ind3)-q_droite(ind3));
	pq4=(p_haut(ind4)-p_bas(ind4)-q(ind4)-q_droite(ind4))/6;
	pq5=(-p(ind5)-p_bas(ind5)+q_gauche(ind5)-q_droite(ind5))/6;
	pq6=0.25*(-p(ind6)-p_bas(ind6)-q(ind6)-q_droite(ind6));
	pq7=0.25*(p_haut(ind7)+p(ind7)+q_gauche(ind7)+q(ind7));
	pq8=0.5*(p_haut(ind8)+p(ind8));
	pq9=0.5*(q_gauche(ind9)+q(ind9));
	
	if(trace)
		h=figure();
		h2=figure();
	end
	
	for it=1:it_max

		u_prec = u;
	
		u_haut=[zeros(1,ncols);u(1:end-1,:)];
		u_bas =[u(2:end,:);zeros(1,ncols)];
		u_gauche=[zeros(nrows,1),u(:,1:end-1)];
		u_droite=[u(:,2:end),zeros(nrows,1)];		
		u_bd=u_bas+u_droite;
		u_hg=u_haut+u_gauche;		
		u3=0.25*(u_bd(ind3)+u_hg(ind3));
		u4=(u_bd(ind4)+u_haut(ind4))/3;
		u5=(u_bd(ind5)+u_gauche(ind5))/3;
		u6=0.5*u_bd(ind6);
		u7=0.5*u_hg(ind7);
		u8=u_haut(ind8);
		u9=u_gauche(ind9);
		
		u(ind3)=u3+pq3;
		u(ind4)=u4+pq4;
		u(ind5)=u5+pq5;
		u(ind6)=u6+pq6;
		u(ind7)=u7+pq7;
		u(ind8)=u8+pq8;
		u(ind9)=u9+pq9;
		
		% Cas particuliers
		u(end,end)=0.5*(u(end-1,end)+u(end,end-1));

		rel_res = norm(u_prec(imask)-u(imask))/norm(u_prec(imask));
		if( rel_res < tol )
			disp('Convergence reached')
			break;
		end
		
		if(save_rmse)
			moyenne_ecarts=mean(u(omega>0)-ground_truth(omega>0));
			u=u-moyenne_ecarts;
			rmse=sqrt(sum((u(omega>0)-ground_truth(omega>0)).^2)/npix);
			tab_rmse=[tab_rmse,rmse];			
		end
		
		if(trace)	
			% Affichage de la surface
			if(mod(it,100)==1)
				disp(sprintf('it %d - rel. res : %.9f',it,rel_res));
				figure(h)
				surfl((u),[0 90])
				axis ij
				view(-45,15)				
				axis equal
				shading flat
				colormap gray
				
				if(save_rmse)
					figure(h2)
					plot((1:it)/size(p,1),tab_rmse)
					xlabel('$k/n$','Interpreter','Latex','FontSize',28)
					ylabel('RMSE','Interpreter','Latex','FontSize',28)
					set(gca,'FontSize',18)
				end
				
			end
		end		
	end
	
	

end
