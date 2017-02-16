function [b,uinf] = femeg_indep_fs(p,t,pos,q,D,varargin)
% FEMEG_INDEP_FS computes the source (independent) vector for the EEG-FP
%                using the full subtraction approach
%
% [b, uinf] = FEMEG_INDEP_FS( p, t, pos, q, D, n) returns the independent
% vector b and the solution in an homogeneous and infinite medium uinf
% provided a mesh with nodes p and t, an electrical conductivity tensor 
% field D, and a source located at pos with moment q. The last entry 
% (optional) refers to the order of numerical integration scheme used (n=7 
% by default) This is done in accordance to: 
%
%   + Beltrachini, L., "A finite element solution of the EEG forward problem 
%     for multipolar sources", Submitted
%
% Inputs:
%    p: nodes defining the tetrahedral mesh (size: Np x 3).
%    t: elements defining the tetrahedral mesh (size: Nt x No+1, with No:
%       number of nodes for a given basis function, i.e. No=4 for first 
%       order basis functions, and No=10 for second order basis functions).
%       The last column is used to group elements belonging to the same
%       layer (by using the same number).
%    pos: Posiotion of the multipolar source (size: 1 x 3).
%    q: Multipolar moment of the source (size: 1 x Nm, with Nm: number of
%       elements defining the multipolar moment; Nm = 1 for monopoles, 
%       Nm = 3 for dipoles, and Nm = 6 for quadrupoles, defined as 
%       [Qxx Qxy Qxz Qyy Qyz Qzz]).
%    D: Electrical conductivity tensor for each element (size: Nt x 6, in
%       the form [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]).
%    n: (optional) Numerical integration order (7 by default).
%
% Outputs:
%    b: Source (independent) vector (size: Np x 1).
%    uinf: Solution of the corresponding source in free space (size: Np x 1).
%
% The function femeg_som can be used for defining a second order mesh as a
% function of a first order mesh.
%
% See also FEMEG_SOM

%   This file is part of the FEMEG toolbox.
%   Author: Leandro Beltrachini <BeltrachiniL at cardiff.ac.uk>


%% General considerations

% Set integration order
if nargin==6, n=varargin{1};else n=7; end

% Set electrical conductivity tensors as needed
Tpos=tsearchn(p,t(:,1:4),pos);
Dy=D(Tpos,1);D(:,[1,4,6])=Dy-D(:,[1,4,6]);D(:,[2,3,5])=-D(:,[2,3,5]);

%% Compute volume integrals 
b1=femeg_int_b1(p,t,q,pos,n,Dy,D);

%% Compute surface integrals
b2=femeg_int_b2(p,t,q,pos,n);

%% Compute final independent vector
b=b1+b2;

%% Compute solution infinite space (isotropic medium)
uinf=p-repmat(pos,size(p,1),1);
nui2=sum(uinf.^2,2);

if numel(q)==1 % monopole
    uinf=q/4/pi/Dy./sqrt(nui2);
elseif numel(q)==3 % dipole
    uinf=1/4/pi/Dy*(q(1)*uinf(:,1)+q(2)*uinf(:,2)+q(3)*uinf(:,3))./(nui2.^(3/2));
elseif numel(q)==6 % quadrupole
    uinf=-1/4/pi/Dy/2*(nui2*sum(q([1,4,6]))-3*(q(1)*uinf(:,1).^2+q(4)*uinf(:,2).^2+...
          q(6)*uinf(:,3).^2+2*q(2)*uinf(:,1).*uinf(:,2)+2*q(3)*uinf(:,1).*uinf(:,3)+...
          2*q(5)*uinf(:,2).*uinf(:,3)))./(nui2.^(5/2));
end

end

%% Compute integrals over tetrahedrons

function [b1]=femeg_int_b1(p,t,q,pos,n,Dy,D)

if size(t,2)<=5, ord=1; va=3; else ord=2; va=10; end

% Get nodes and weights
[t0,w0]=jacpts(n,0,0,[0,1]);
[t1,w1]=jacpts(n,1,0,[0,1]);
[t2,w2]=jacpts(n,2,0,[0,1]);
cont=1;rn=zeros(3,n^3);wt=zeros(n^3,1);
for ii=1:n
    for jj=1:n
        for kk=1:n
            rn(:,cont)=[t2(kk);t1(jj)*(1-t2(kk));t0(ii)*(1-t1(jj))*(1-t2(kk))];
            wt(cont)=w2(kk)*w1(jj)*w0(ii);
            cont=cont+1;
        end
    end
end

% Find elements belonging to regions with different conductivity than Dy
Ind_cn=sum(abs(D),2)>1e-12; % Ind_cn=1:size(D,1);

% Copute volume coordintes
[b,c,d,V,a]=femeg_vol_coord(p,t(Ind_cn,1:4));

% Compute integral
gu=zeros(nnz(Ind_cn),va);
for inn=1:n^3
    
    % Transformed points
    xn=rn(1,inn)*(p(t(Ind_cn,1),:)-p(t(Ind_cn,4),:))+...
       rn(2,inn)*(p(t(Ind_cn,2),:)-p(t(Ind_cn,4),:))+...
       rn(3,inn)*(p(t(Ind_cn,3),:)-p(t(Ind_cn,4),:))+...
       p(t(Ind_cn,4),:);
    x=xn-repmat(pos,size(xn,1),1);nx2=sum(x.^2,2);
   
    if numel(q)==1 % Monopolar source
        In=-q*x./repmat(nx2.^(3/2),1,3);
        
    elseif numel(q)==3 % Dipolar source
        In=[(q(1)*nx2-3*(q(1)*x(:,1)+q(2)*x(:,2)+q(3)*x(:,3)).*x(:,1))./(nx2.^(5/2)),...
            (q(2)*nx2-3*(q(1)*x(:,1)+q(2)*x(:,2)+q(3)*x(:,3)).*x(:,2))./(nx2.^(5/2)),...
            (q(3)*nx2-3*(q(1)*x(:,1)+q(2)*x(:,2)+q(3)*x(:,3)).*x(:,3))./(nx2.^(5/2))];
         
    elseif numel(q)==6 % Quadrupolar source
        gu_prov_x=((2*x(:,1)*sum(q([1,4,6]))-3*(q(1)*2*x(:,1)+q(2)*2*x(:,2)+q(3)*2*x(:,3))).*nx2-...
                   (sum(q([1,4,6])).*nx2-3*(q(1)*x(:,1).^2+q(4)*x(:,2).^2+q(6)*x(:,3).^2+2*q(2)*x(:,1).*x(:,2)+...
                   2*q(3)*x(:,1).*x(:,3)+2*q(5)*x(:,2).*x(:,3)))*5.*x(:,1))./(nx2.^(7/2));
        gu_prov_y=((2*x(:,2)*sum(q([1,4,6]))-3*(q(2)*2*x(:,1)+q(4)*2*x(:,2)+q(5)*2*x(:,3))).*nx2-...
                   (sum(q([1,4,6])).*nx2-3*(q(1)*x(:,1).^2+q(4)*x(:,2).^2+q(6)*x(:,3).^2+2*q(2)*x(:,1).*x(:,2)+...
                   2*q(3)*x(:,1).*x(:,3)+2*q(5)*x(:,2).*x(:,3)))*5.*x(:,2))./(nx2.^(7/2));
        gu_prov_z=((2*x(:,3)*sum(q([1,4,6]))-3*(q(3)*2*x(:,1)+q(5)*2*x(:,2)+q(6)*2*x(:,3))).*nx2-...
                   (sum(q([1,4,6])).*nx2-3*(q(1)*x(:,1).^2+q(4)*x(:,2).^2+q(6)*x(:,3).^2+2*q(2)*x(:,1).*x(:,2)+...
                   2*q(3)*x(:,1).*x(:,3)+2*q(5)*x(:,2).*x(:,3)))*5.*x(:,3))./(nx2.^(7/2));       
        In=-[gu_prov_x,gu_prov_y,gu_prov_z]/2;
    end
    
    % Shape functions
    if ord==1
        gu=gu+In.*wt(inn); 

    elseif ord==2
        
        % Compute bbb=gu*D*Lambda
        bbb=[In(:,1).*(b(:,1).*D(Ind_cn,1) + c(:,1).*D(Ind_cn,2) + d(:,1).*D(Ind_cn,3)) + In(:,2).*(b(:,1).*D(Ind_cn,2) + c(:,1).*D(Ind_cn,4) + d(:,1).*D(Ind_cn,5)) + In(:,3).*(b(:,1).*D(Ind_cn,3) + c(:,1).*D(Ind_cn,5) + d(:,1).*D(Ind_cn,6)),...
             In(:,1).*(b(:,2).*D(Ind_cn,1) + c(:,2).*D(Ind_cn,2) + d(:,2).*D(Ind_cn,3)) + In(:,2).*(b(:,2).*D(Ind_cn,2) + c(:,2).*D(Ind_cn,4) + d(:,2).*D(Ind_cn,5)) + In(:,3).*(b(:,2).*D(Ind_cn,3) + c(:,2).*D(Ind_cn,5) + d(:,2).*D(Ind_cn,6)),...
             In(:,1).*(b(:,3).*D(Ind_cn,1) + c(:,3).*D(Ind_cn,2) + d(:,3).*D(Ind_cn,3)) + In(:,2).*(b(:,3).*D(Ind_cn,2) + c(:,3).*D(Ind_cn,4) + d(:,3).*D(Ind_cn,5)) + In(:,3).*(b(:,3).*D(Ind_cn,3) + c(:,3).*D(Ind_cn,5) + d(:,3).*D(Ind_cn,6)),...
             In(:,1).*(b(:,4).*D(Ind_cn,1) + c(:,4).*D(Ind_cn,2) + d(:,4).*D(Ind_cn,3)) + In(:,2).*(b(:,4).*D(Ind_cn,2) + c(:,4).*D(Ind_cn,4) + d(:,4).*D(Ind_cn,5)) + In(:,3).*(b(:,4).*D(Ind_cn,3) + c(:,4).*D(Ind_cn,5) + d(:,4).*D(Ind_cn,6))];
        
        % Volume coordinates
        xi=repmat(xn(:,1),1,4).*b+repmat(xn(:,2),1,4).*c+repmat(xn(:,3),1,4).*d+a;
        xi=xi./repmat(6*V,1,4);
     
        % Compute In=bbb*gxi
        In=[bbb(:,1).*(4*xi(:,1)-1),bbb(:,2).*(4*xi(:,2)-1),bbb(:,3).*(4*xi(:,3)-1),bbb(:,4).*(4*xi(:,4)-1),...
            4*(xi(:,2).*bbb(:,1)+xi(:,1).*bbb(:,2)),4*(xi(:,3).*bbb(:,1)+xi(:,1).*bbb(:,3)),4*(xi(:,4).*bbb(:,1)+xi(:,1).*bbb(:,4)),...
            4*(xi(:,2).*bbb(:,3)+xi(:,3).*bbb(:,2)),4*(xi(:,2).*bbb(:,4)+xi(:,4).*bbb(:,2)),4*(xi(:,4).*bbb(:,3)+xi(:,3).*bbb(:,4))];
        
        gu=gu+In.*wt(inn); 
    end    

end

% explicar
dR=abs((p(t(Ind_cn,1),1)-p(t(Ind_cn,4),1)).*(p(t(Ind_cn,2),2)-p(t(Ind_cn,4),2)).*(p(t(Ind_cn,3),3)-p(t(Ind_cn,4),3))+...
       (p(t(Ind_cn,1),2)-p(t(Ind_cn,4),2)).*(p(t(Ind_cn,2),3)-p(t(Ind_cn,4),3)).*(p(t(Ind_cn,3),1)-p(t(Ind_cn,4),1))+...
       (p(t(Ind_cn,1),3)-p(t(Ind_cn,4),3)).*(p(t(Ind_cn,2),1)-p(t(Ind_cn,4),1)).*(p(t(Ind_cn,3),2)-p(t(Ind_cn,4),2))-...
       (p(t(Ind_cn,3),1)-p(t(Ind_cn,4),1)).*(p(t(Ind_cn,2),2)-p(t(Ind_cn,4),2)).*(p(t(Ind_cn,1),3)-p(t(Ind_cn,4),3))-...
       (p(t(Ind_cn,3),2)-p(t(Ind_cn,4),2)).*(p(t(Ind_cn,2),3)-p(t(Ind_cn,4),3)).*(p(t(Ind_cn,1),1)-p(t(Ind_cn,4),1))-...
       (p(t(Ind_cn,2),1)-p(t(Ind_cn,4),1)).*(p(t(Ind_cn,1),2)-p(t(Ind_cn,4),2)).*(p(t(Ind_cn,3),3)-p(t(Ind_cn,4),3)));

% 
be1=gu.*repmat(dR,1,size(gu,2))/4/pi/Dy;

if size(t,2)<=5
    be1=[be1(:,1).*(b(:,1).*D(Ind_cn,1) + c(:,1).*D(Ind_cn,2) + d(:,1).*D(Ind_cn,3)) + be1(:,2).*(b(:,1).*D(Ind_cn,2) + c(:,1).*D(Ind_cn,4) + d(:,1).*D(Ind_cn,5)) + be1(:,3).*(b(:,1).*D(Ind_cn,3) + c(:,1).*D(Ind_cn,5) + d(:,1).*D(Ind_cn,6)),...
         be1(:,1).*(b(:,2).*D(Ind_cn,1) + c(:,2).*D(Ind_cn,2) + d(:,2).*D(Ind_cn,3)) + be1(:,2).*(b(:,2).*D(Ind_cn,2) + c(:,2).*D(Ind_cn,4) + d(:,2).*D(Ind_cn,5)) + be1(:,3).*(b(:,2).*D(Ind_cn,3) + c(:,2).*D(Ind_cn,5) + d(:,2).*D(Ind_cn,6)),...
         be1(:,1).*(b(:,3).*D(Ind_cn,1) + c(:,3).*D(Ind_cn,2) + d(:,3).*D(Ind_cn,3)) + be1(:,2).*(b(:,3).*D(Ind_cn,2) + c(:,3).*D(Ind_cn,4) + d(:,3).*D(Ind_cn,5)) + be1(:,3).*(b(:,3).*D(Ind_cn,3) + c(:,3).*D(Ind_cn,5) + d(:,3).*D(Ind_cn,6)),...
         be1(:,1).*(b(:,4).*D(Ind_cn,1) + c(:,4).*D(Ind_cn,2) + d(:,4).*D(Ind_cn,3)) + be1(:,2).*(b(:,4).*D(Ind_cn,2) + c(:,4).*D(Ind_cn,4) + d(:,4).*D(Ind_cn,5)) + be1(:,3).*(b(:,4).*D(Ind_cn,3) + c(:,4).*D(Ind_cn,5) + d(:,4).*D(Ind_cn,6))];
end

% Assemble vector
b1=zeros(size(p,1),1);
for kk=1:size(be1,2)
    b1=b1+accumarray(t(Ind_cn,kk),be1(:,kk)./V/6,[size(p,1),1]);
end

end


%% Integration over triangles

function [b2,tri]=femeg_int_b2(p,t,q,pos,n)

% Get nodes and weights
[tj,wj]=jacpts(n,1,0,[0,1]);
[tl,wl]=legpts(n,[0,1]);
cont=1;rn=zeros(2,n^2);wt=zeros(n^2,1);
for ii=1:n
    for jj=1:n
        rn(:,cont)=[tj(ii);(1-tj(ii))*tl(jj)];
        wt(cont)=wj(ii)*wl(jj);
        cont=cont+1;
    end
end

% Find surface triangles
tri=surftri(p,t);

% Compute Namt=dx(s,t)/ds x dx(s,t)/dt
Namt=cross(p(tri(:,2),:)-p(tri(:,1),:),p(tri(:,3),:)-p(tri(:,1),:),2);

% Compute integral
gu=zeros(size(tri,1),size(tri,2));
for inn=1:n^2
    
    % Transformed points
    x=rn(1,inn)*(p(tri(:,2),:)-p(tri(:,1),:))+rn(2,inn)*(p(tri(:,3),:)-p(tri(:,1),:))+p(tri(:,1),:);
    x=x-repmat(pos,size(x,1),1);nx2=sum(x.^2,2);
   
    % Compute gu_prov=gu*Namt in integration nodes
    if numel(q)==1 % For a monopolar source
        gu_prov=-q*dot(x,Namt,2)./(nx2.^(3/2));
    elseif numel(q)==3 % For a dipolar source
        gu_prov=((q(1)*Namt(:,1)+q(2)*Namt(:,2)+q(3)*Namt(:,3)).*nx2-3*(q(1)*x(:,1)+q(2)*x(:,2)+q(3)*x(:,3)).*dot(x,Namt,2))./(nx2.^(5/2));
    elseif numel(q)==6 % For a quadrupolar source
        gu_prov_x=((2*x(:,1)*sum(q([1,4,6]))-3*(q(1)*2*x(:,1)+q(2)*2*x(:,2)+q(3)*2*x(:,3))).*nx2-...
                   (sum(q([1,4,6])).*nx2-3*(q(1)*x(:,1).^2+q(4)*x(:,2).^2+q(6)*x(:,3).^2+2*q(2)*x(:,1).*x(:,2)+...
                   2*q(3)*x(:,1).*x(:,3)+2*q(5)*x(:,2).*x(:,3)))*5.*x(:,1))./(nx2.^(7/2));
        gu_prov_y=((2*x(:,2)*sum(q([1,4,6]))-3*(q(2)*2*x(:,1)+q(4)*2*x(:,2)+q(5)*2*x(:,3))).*nx2-...
                   (sum(q([1,4,6])).*nx2-3*(q(1)*x(:,1).^2+q(4)*x(:,2).^2+q(6)*x(:,3).^2+2*q(2)*x(:,1).*x(:,2)+...
                   2*q(3)*x(:,1).*x(:,3)+2*q(5)*x(:,2).*x(:,3)))*5.*x(:,2))./(nx2.^(7/2));
        gu_prov_z=((2*x(:,3)*sum(q([1,4,6]))-3*(q(3)*2*x(:,1)+q(5)*2*x(:,2)+q(6)*2*x(:,3))).*nx2-...
                   (sum(q([1,4,6])).*nx2-3*(q(1)*x(:,1).^2+q(4)*x(:,2).^2+q(6)*x(:,3).^2+2*q(2)*x(:,1).*x(:,2)+...
                   2*q(3)*x(:,1).*x(:,3)+2*q(5)*x(:,2).*x(:,3)))*5.*x(:,3))./(nx2.^(7/2));       
        gu_prov=dot(-[gu_prov_x,gu_prov_y,gu_prov_z],Namt,2)/2;
    end
    
    % Shape functions
    if size(tri,2)==3
        gu_prov=kron([(1-rn(1,inn)-rn(2,inn)),rn(1,inn),rn(2,inn)],gu_prov);
    elseif size(tri,2)==6
        gu_prov=kron([(1-rn(1,inn)-rn(2,inn)).*(1-2*rn(1,inn)-2*rn(2,inn)),...
                 rn(1,inn).*(2*rn(1,inn)-1),...
                 rn(2,inn).*(2*rn(2,inn)-1),...
                 4*rn(2,inn).*(1-rn(1,inn)-rn(2,inn)),...
                 4*rn(1,inn).*(1-rn(1,inn)-rn(2,inn)),...
                 4*rn(1,inn).*rn(2,inn)],gu_prov);
    end
            
    % Final value 
    gu=gu+gu_prov.*wt(inn);
end

% Assemble vector
be2=gu/4/pi; 
b2=zeros(size(p,1),1);
for kk=1:size(be2,2)
    b2=b2-accumarray(tri(:,kk),be2(:,kk),[size(p,1),1]);
end
  
end

