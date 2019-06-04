%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A 99 LINE TOPOLOGY OPTIMIZATION CODE BY OLE SIGMUND, OCTOBER 1999
% MODIFIED FOR 3D MULTISCALE DESIGN VIA SURROGATE MODEL, LLNL, JULY 2018
%
% This work was produced under the auspices of the U.S. Department of Energy by
% Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344.
%
% This work was prepared as an account of work sponsored by an agency of the
% United States Government. Neither the United States Government nor Lawrence
% Livermore National Security, LLC, nor any of their employees makes any warranty,
% expressed or implied, or assumes any legal liability or responsibility for the
% accuracy, completeness, or usefulness of any information, apparatus, product, or
% process disclosed, or represents that its use would not infringe privately owned
% rights. Reference herein to any specific commercial product, process, or service
% by trade name, trademark, manufacturer, or otherwise does not necessarily
% constitute or imply its endorsement, recommendation, or favoring by the United
% States Government or Lawrence Livermore National Security, LLC. The views and
% opinions of authors expressed herein do not necessarily state or reflect those
% of the United States Government or Lawrence Livermore National Security, LLC,
% and shall not be used for advertising or product endorsement purposes.
%
% LLNL-CODE-757968
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function top(nelx, nely, nelz, volfrac, rmin, truss, Es, vs, minVF, maxVF, maxit)
% INITIALIZE
x(1:nelx, 1:nely, 1:nelz) = volfrac;
Gs = Es / (2*(1+vs));
loop = 0; change = 1.0;
nnx = nelx+1; nny = nely+1; nnz = nelz+1;
colormap(gray); caxis([0.0, 1.0]);
% START ITERATION
while change > 0.01 && loop < maxit
	loop = loop + 1;
	xold = x;
	% FE-ANALYSIS
	[U] = FE(nelx, nely, nelz, nnx, nny, nnz, x, truss, Es, vs, Gs);
	% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
	c = 0.0; dc = zeros(nelx, nely, nelz);
	for elz = 1:nelz; for ely = 1:nely for elx = 1:nelx
		[KE]   = get_KE(truss, x, Es, vs, Gs, elx, ely, elz, 0);
		[DKE]  = get_KE(truss, x, Es, vs, Gs, elx, ely, elz, 1);
		[dofs] = get_elem_dofs(nnx, nny, nnz, elx, ely, elz);
		Ue = U(dofs,1);
		c = c + Ue'*KE*Ue;
		dc(elx,ely,elz) = -Ue'*DKE*Ue;
	end; end; end
	% FILTERING OF SENSITIVITIES
	[dc]   = check(nelx, nely, nelz, rmin, x, dc);
	% DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
	[x]    = OC(nelx, nely, nelz, x, volfrac, dc, minVF, maxVF);
	% PRINT RESULTS
	change = max(max(max(abs(x-xold))));
	disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',c) ...
       ' Vol.: ' sprintf('%6.3f',sum(sum(sum(x)))/(nelx*nely*nelz)) ...
        ' ch.: ' sprintf('%6.3f',change )])
	% PLOT DENSITIES
	viz3d(nelx, nely, nelz, x, volfrac, nelx==1);
	% SAVE PARAMETER VALUES (ELEMENT DENSITIES AND ROD DIAMETERS)
	xOut = reshape(x,[],1);             save('-ascii','elVolFrac.txt', 'xOut');
	dOut = reshape(get_d(truss,x),[],1); save('-ascii','elRodDiam.txt','dOut');
end
%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew] = OC(nelx, nely, nelz, x, volfrac, dc, minVF, maxVF)
l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
	lmid = 0.5*(l2 + l1);
	xnew = max(minVF, max(x-move, min(maxVF, min(x+move,x.*sqrt(-dc./lmid)))));
	if sum(sum(sum(xnew))) - volfrac*nelx*nely*nelz > 0;
		l1 = lmid;
	else
		l2 = lmid;
	end
end
%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn] = check(nelx, nely, nelz, rmin, x, dc)
dcn=zeros(size(dc));
for elz = 1:nelz; for ely = 1:nely; for elx = 1:nelx
	sum = 0.0;
	for k = max(elz-round(rmin),1):min(elz+round(rmin),nelz)
		for j = max(ely-round(rmin),1):min(ely+round(rmin),nely)
			for i = max(elx-round(rmin),1):min(elx+round(rmin),nelx)
				fac = rmin - sqrt((elx-i)^2+(ely-j)^2+(elz-k)^2);
				sum = sum + max(0,fac);
				dcn(elx,ely,elz) = dcn(elx,ely,elz) + max(0,fac)*x(i,j,k)*dc(i,j,k);
			end
		end
	end
	dcn(elx,ely,elz) = dcn(elx,ely,elz) / (x(elx,ely,elz)*sum);
end; end; end
%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [U] = FE(nelx, nely, nelz, nnx, nny, nnz, x, truss, Es, vs, Gs)
K = sparse(3*nnx*nny*nnz, 3*nnx*nny*nnz);
F = sparse(3*nnx*nny*nnz,1);
U = sparse(3*nnx*nny*nnz,1);
for elz = 1:nelz; for ely = 1:nely; for elx = 1:nelx
	[KE]   = get_KE(truss, x, Es, vs, Gs, elx, ely, elz, 0);
	[dofs] = get_elem_dofs(nnx, nny, nnz, elx, ely, elz);
	K(dofs,dofs) = K(dofs,dofs) + KE;
end; end; end
% DEFINE LOADS AND SUPPORTS(HALF MBB-BEAM)
coords = zeros(nnx*nny*nnz,3);
n = 0;
for k = 1:nnz; for j = 1:nny; for i = 1:nnx
	n = n+1; coords(n,1) = i-1; coords(n,2) = j-1; coords(n,3) = k-1;
end; end; end
midplane_nodes = find(coords(:,2)==0);
loaded_nodes   = intersect(find(coords(:,3)==nelz), find(coords(:,2)==0));
fixed_nodes    = intersect(find(coords(:,3)==0), find(coords(:,2)==nely));
fixeddofs      = zeros(size(midplane_nodes,1) + 2*size(fixed_nodes,1),1);
for i = loaded_nodes'; F(3*(i-1)+3) = -1.0/nnx; end
n = 1;
for i = midplane_nodes'; for j=[2]; fixeddofs(n,1) = 3*(i-1)+j; n =n+1; end; end
for i = fixed_nodes'; for j=[1,3];  fixeddofs(n,1) = 3*(i-1)+j; n =n+1; end; end
alldofs   = [1:3*nnx*nny*nnz];
freedofs  = setdiff(alldofs,fixeddofs);
% SOLVING
U(freedofs,:)  = K(freedofs,freedofs) \ F(freedofs,1);
U(fixeddofs,:) = 0;
%%%%%%%%% ELEMENT AND NODE NUMBERING IN 3D MESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [num] = get_num(nx, ny, nz, i, j, k)
num = (nx*ny)*(k-1) + nx*(j-1) + i;
%%%%%%%%% GLOBAL DOFS FOR A GIVEN ELEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dofs] = get_elem_dofs(nnx, nny, nnz, elx, ely, elz)
n = get_num(nnx, nny, nnz, elx, ely, elz);
N = [n; n+1; n+nnx+1; n+nnx; n+nnx*nny; n+nnx*nny+1;
     n+nnx*nny+nnx+1; n+nnx*nny+nnx];
dofs = zeros(24,1); 
for j = 1:8; for i = 1:3; dofs(3*(j-1)+i) = 3*(N(j)-1)+i; end; end;
%%%%%%%%% INTEGRATE ELASTICITY TENSOR CE TO GET KE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE] = get_KE(truss, x, Es, vs, Gs, i, j, k, deriv)
KE = zeros(24,24);
CE = get_CE(truss, x, Es, vs, Gs, i, j, k, deriv);
for l = 1:8
	r = (sqrt(3)/3) * (-1 + 2*any([2,3,6,7]==l)); rp = (1+r); rm = (1-r);
	s = (sqrt(3)/3) * (-1 + 2*any([3,4,7,8]==l)); sp = (1+s); sm = (1-s);
	t = (sqrt(3)/3) * (-1 + 2*any([5,6,7,8]==l)); tp = (1+t); tm = (1-t);
	DN = [-sm*tm, -rm*tm, -rm*sm;    sm*tm, -rp*tm, -rp*sm;
           sp*tm,  rp*tm, -rp*sp;   -sp*tm,  rm*tm, -rm*sp;
          -sm*tp, -rm*tp,  rm*sm;    sm*tp, -rp*tp,  rp*sm;
           sp*tp,  rp*tp,  rp*sp;   -sp*tp, rm*tp,  rm*sp] / 8;
	B = DN * 2*eye(3); G = kron(B', eye(3)); KE = KE + G' * CE * G / 4;
end
%%%%%%%%% DEFINE ELASTICITY TENSOR FOR DIFFERENT TRUSSES %%%%%%%%%%%%%%%%%%%%%%%
function [CE] = get_CE(truss, x, Es, vs, Gs, i, j, k, D)
p = x(i, j, k);
if strcmpi(truss, 'iso');       TM = @(p,Es,vs,Gs,D) iso_moduli(p,Es,vs,Gs,D);
elseif strcmpi(truss, 'octet'); TM = @(p,Es,vs,Gs,D) octet_moduli(p,Es,vs,Gs,D);
elseif strcmpi(truss, 'orc');   TM = @(p,Es,vs,Gs,D) orc_moduli(p,Es,vs,Gs,D);
elseif strcmpi(truss, 'bound'); TM = @(p,Es,vs,Gs,D) bound_moduli(p,Es,vs,Gs,D);
else;                           TM = @(p,Es,vs,Gs,D) simp_moduli(p,Es,vs,Gs,D);
end
[E, v, G] = TM(p,Es,vs,Gs,0);  if D; [DE, Dv, DG] = TM(p,Es,vs,Gs,1); end
if D == 0
	C1111 = E * (1.0 - v) / (1.0 - v - 2*v^2);
	C1122 = (E * v) / (1.0 - v - 2*v^2);
	C1212 = G;
else % return the deriviatives instead
	C1111 = ((DE*(1-v)-E*Dv)*(1-v-2*v^2)-E*(1-v)*(-Dv-4*v*Dv)) / (1-v-2*v^2)^2;
	C1122 = ((DE*v+E*Dv)*(1-v-2*v^2)-E*(1-v)*(-Dv-4*v*Dv)) / (1-v-2*v^2)^2;
	C1212 = DG;
end
CE = [C1111   0     0     0   C1122   0     0     0   C1122;
        0   C1212   0   C1212   0     0     0     0     0  ;
        0     0   C1212   0     0     0   C1212   0     0  ;
        0   C1212   0   C1212   0     0     0     0     0  ;
      C1122   0     0     0   C1111   0     0     0   C1122;
        0     0     0     0     0   C1212   0   C1212   0  ;
        0     0   C1212   0     0     0   C1212   0     0  ;
        0     0     0     0     0   C1212   0   C1212   0  ;
      C1122   0     0     0   C1122   0     0     0   C1111];
%%%%%%%%% TRUSS-SPECIFIC MECHANICS MODELS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E,v,G] = iso_moduli(p, Es, vs, Gs, deriv)
E = Es * ( 2.05292e-01 - 3.30265e-02*vs) * (p^(1-deriv)) * (1+0*deriv) + ...
	     ( 8.12145e-02 + 2.72431e-01*vs) * (p^(2-deriv)) * (1+1*deriv) + ...
	     ( 6.49737e-01 - 2.42374e-01*vs) * (p^(3-deriv)) * (1+2*deriv);
v =      ( 2.47760e-01 + 1.69804e-02*vs) * (1-deriv) + ...
	     (-1.59293e-01 + 7.38598e-01*vs) * (p^(1-deriv)) * (1+0*deriv) + ...
	     (-1.86279e-01 - 4.83229e-01*vs) * (p^(2-deriv)) * (1+1*deriv) + ...
	     ( 9.77457e-02 + 7.26595e-01*vs) * (p^(3-deriv)) * (1+2*deriv);
G = Gs * ( 1.63200e-01 + 1.27910e-01*vs) * (p^(1-deriv)) * (1+0*deriv) + ...
	     ( 6.00810e-03 + 4.13331e-01*vs) * (p^(2-deriv)) * (1+1*deriv) + ...
	     ( 7.22847e-01 - 3.56032e-01*vs) * (p^(3-deriv)) * (1+2*deriv);
function [E,v,G] = octet_moduli(p, Es, vs, Gs, deriv)
E = Es * ( 1.36265e-01 - 1.22204e-02*vs) * (p^(1-deriv)) * (1+0*deriv) + ...
	     ( 8.57991e-02 + 6.63677e-02*vs) * (p^(2-deriv)) * (1+1*deriv) + ...
	     ( 7.39887e-01 - 6.26129e-02*vs) * (p^(3-deriv)) * (1+2*deriv);
v =      ( 3.29529e-01 + 1.86038e-02*vs) * (1-deriv) + ...
	     (-1.42155e-01 + 4.57806e-01*vs) * (p^(1-deriv)) * (1+0*deriv) + ...
	     (-3.29837e-01 + 5.59823e-02*vs) * (p^(2-deriv)) * (1+1*deriv) + ...
	     ( 1.41233e-01 + 4.72695e-01*vs) * (p^(3-deriv)) * (1+2*deriv);
G = Gs * ( 2.17676e-01 + 7.22515e-02*vs) * (p^(1-deriv)) * (1+0*deriv) + ...
	     (-7.63847e-02 + 1.31601e+00*vs) * (p^(2-deriv)) * (1+1*deriv) + ...
	     ( 9.11800e-01 - 1.55261e+00*vs) * (p^(3-deriv)) * (1+2*deriv);
function [E,v,G] = orc_moduli(p, Es, vs, Gs, deriv)
E = Es * ( 1.34332e-01 - 7.06384e-02*vs) * (p^(1-deriv)) * (1+0*deriv) + ...
	     ( 2.59957e-01 + 8.51515e-01*vs) * (p^(2-deriv)) * (1+1*deriv) + ...
	     ( 6.53902e-01 - 7.29803e-01*vs) * (p^(3-deriv)) * (1+2*deriv);
v =      ( 3.38525e-01 + 7.04361e-03*vs) * (1-deriv) + ...
	     (-4.25721e-01 + 4.14882e-01*vs) * (p^(1-deriv)) * (1+0*deriv) + ...
	     (-7.68215e-02 + 5.58948e-01*vs) * (p^(2-deriv)) * (1+1*deriv) + ...
	     ( 1.64073e-01 + 3.98374e-02*vs) * (p^(3-deriv)) * (1+2*deriv);
G = Gs * ( 1.96762e-01 + 1.66705e-01*vs) * (p^(1-deriv)) * (1+0*deriv) + ...
	     ( 1.30938e-01 + 1.72565e-01*vs) * (p^(2-deriv)) * (1+1*deriv) + ...
	     ( 6.45455e-01 - 2.87424e-01*vs) * (p^(3-deriv)) * (1+2*deriv);
function [E,v,G] = bound_moduli(p, Es, vs, Gs, deriv)
Ks = 1.0 / (3*(1-2*vs));
K = Ks + (1-p) / ( -1.0/Ks + p/(Ks + (4.0*Gs)/3.0) );
G = Gs + (1-p) / ( -1.0/Gs + (2.0*p*(Ks+2.0*Gs)) / (5.0*Gs*(Ks+(4.0*Gs)/3.0)) );
E = 9*K*G/(3*K+G);
v = (3*K-2*G) / (2*(3*K+G));
if deriv
	DK = (p - 1)/(((4*Gs)/3 + Ks)*(p/((4*Gs)/3 + Ks) - 1/Ks)^2) - ...
		1/(p/((4*Gs)/3 + Ks) - 1/Ks);
	DG = 1/(1/Gs - (2*p*(2*Gs + Ks))/(5*Gs*((4*Gs)/3 + Ks))) + ...
		(2*(2*Gs + Ks)*(p - 1))/(5*Gs*((4*Gs)/3 + Ks)*(1/Gs - ...
		(2*p*(2*Gs + Ks))/(5*Gs*((4*Gs)/3 + Ks)))^2);
	DE = ( 9*(3*K+G)*(DK*G+K*DG) - 9*K*G*(3*DK+DG) ) / (3*K+G)^2;
	Dv = ( 2*(3*K+G)*(3*DK-2*DG) - 2*(3*K-2*G)*(3*DK+DG) ) / (2*(3*K+G))^2;
	G = DG;
	E = DE;
	v = Dv;
end
function [E,v,G] = simp_moduli(p, Es, vs, Gs, deriv)
E = Es * p^(3-deriv) * (1+2*deriv);
v = vs * (1-deriv);
G = Gs * p^(3-deriv) * (1+2*deriv);
%%%%%%%%% TRUSS-SPECIFIC ROD DIAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [d] = get_d(truss, p)
if strcmpi(truss, 'iso')
	 d = 2.04920e-02 + 1.05076e+00*p - 1.59468e+00*(p.^2) + 1.09799e+00*(p.^3);
elseif strcmpi(truss, 'octet')
	d = 1.64505e-02 + 9.23773e-01*p - 1.61345e+00*(p.^2) + 1.23729e+00*(p.^3);
elseif strcmpi(truss, 'orc')
	d = 2.32950e-02 + 1.31602e+00*p - 2.28842e+00*(p.^2) + 1.90225e+00*(p.^3);
else
	 d = -1*ones(size(p));
end
%%%%%%%%% 3D VISUALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function viz3d(nelx, nely, nelz, x, volfrac, is2D)
y = zeros(nelx+2, nely+2, nelz+2); y(2:nelx+1, 2:nely+1, 2:nelz+1) = x;
if is2D; T=0; A=90; E=0; else; T=volfrac; A=142.5; E=30; end;
nf = nelx*nely*(nelz+1) + nelx*(nely+1)*nelz + (nelx+1)*nely*nelz; n = 0;
X = zeros(4,nf); Y = zeros(4,nf); Z = zeros(4,nf); C = zeros(1,nf);
for k = 1:nelz+1; for j = 1:nely+1; for i = 1:nelx+1;
	I = i-1; J = j-1; K = k-1; L = i+1; M = j+1; N = k+1;
	cz = max(y(L,M,k:N)); cy = max(y(L,j:M,N)); cx = max(y(i:L,M,N));
	dz = min(y(L,M,k:N)); dy = min(y(L,j:M,N)); dx = min(y(i:L,M,N));
	if cz > T && dz < T+is2D; n = n+1; C(1,n) = 1-cz; 
		X(:,n) = [I,i,i,I]'; Y(:,n) = [J,J,j,j]'; Z(:,n) = [K,K,K,K]';
	end
	if cy > T && dy < T+is2D; n = n+1; C(1,n) = 1-cy;
		X(:,n) = [I,i,i,I]'; Y(:,n) = [J,J,J,J]'; Z(:,n) = [K,K,k,k]';
	end
	if cx > T && dx < T+is2D; n = n+1; C(1,n) = 1-cx;
		X(:,n) = [I,I,I,I]'; Y(:,n) = [J,j,j,J]'; Z(:,n) = [K,K,k,k]';
	end
end; end; end
patch(X(:,1:n), Y(:,1:n), Z(:,1:n), C(1,1:n), 'EdgeColor', 'none');
view(A,E); axis equal; axis tight; axis off; pause(1e-3);  
