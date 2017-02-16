function [u] = femeg_anisphere( X, pos, mom, radii, cond, nmax)
%FEMEG_ANISPHERE   Computes the analytical solution of the EEG-FP in any point of a
%                  multilayered spherical conductor model
%
%   u = FEMEG_ANISPHERE(X, pos, mom, radii, cond) returns the electric potential in X 
%   due to a source placed in pos with moment mom considering a multilayered spherical 
%   model with specified radii and conductivities (cond). The solution is based on the 
%   following papers:
%
%
%   + Beltrachini, L., "A finite element solution of the EEG forward problem for 
%      multipolar sources", Submitted
% 
%   + de Munck, J.C., and Peters, M.J. "A fast method to compute the potential in 
%      the multisphere model". IEEE Trans Biomed Eng., 40(11):1166–1174, 1993.
%
%
%   u = FEMEG_ANISPHERE(X, pos, mom, radii, cond, nmax) returns the electric potential
%   considering nmax terms in the series expansion (1000 by default)
%
%
%   Inputs:
%      X: Nodes where to compute the analytical solution (size: Np x 3, with 
%         Np: number of points).
%      pos: Posiotion of the multipolar sources (size: Ns x 3, with 
%           NS: number of sources).
%      mom: Multipolar moment of the sources (size: Ns x Nm, with 
%           Nm: number of elements defining the multipolar moment; Nm = 1 for
%           monopoles, Nm = 3 for dipoles, and Nm = 6 for quadrupoles, defined as 
%           [Qxx Qxy Qxz Qyy Qyz Qzz]).
%      radii: Radii of the spherical layers, from the outermost to the innermost
%             (size: 1 x Nl, with Nl: number of layers).
%      cond: Electrical conductivity of each layer (size: 2 x Nl). The first row 
%            indicates the radial conductivity, while the second indicates the 
%            tangential conductivity (both defined from the outermost to the 
%            innermost layer).
%      nmax: Number of terms considered in the spherical harmonics' series.
%
%
%   Outputs:
%      u: Electric potential for each source (size: Np x Ns).
%
%
%Example: 
%
%%% Head Model (WM,GM,CSF,SKULL,SCALP)
% radii = [0.1,0.095,0.091,0.088,0.063]; % [m]
% cond = [.33,0.001,1.79,.33,.14; ... % radial conductivity [S/m]
%         .33,0.01,1.79,.33,.14];     % tangential conductivity [S/m]
%%% Source
% pos = [0,0,0.07]; % [m]
% mom = [1,0,0]*1e-7; % Dipolar source [A.m]
%
%%% Nodes for computing the potential
% [p,t] = femeg_sphere(0.1,3); % [m]
%
%%% Compute potential
% u = femeg_anisphere(p,pos,mom,radii,cond);
%
%%% Plot results
% figure,vis(p,t,u),colorbar

%   This file is part of the FEMEG toolbox.
%   Author: Nicolas von Ellenrieder and Leandro Beltrachini <BeltrachiniL at cardiff.ac.uk>


if nargin==5, nmax=1000; end
s=[radii;cond];di=[pos,mom];
if size(s,1)==2, s(3,:)=s(2,:); end
nd=size(di,1);
ns=size(X,1);
re=sqrt(X(:,1).^2+X(:,2).^2+X(:,3).^2);
r=[s(1,:),0];
N=size(s,2);
u11=ones(1,nmax);u12=zeros(1,nmax);u21=u12;u22=u11;
for ii=1:N
    [t11,t12,t21,t22]=multiUani(r(ii),r(ii+1),s(2:3,ii),nmax);
    a11=u11.*t11+u12.*t21;
    a12=u11.*t12+u12.*t22;
    a21=u21.*t11+u22.*t21;
    a22=u21.*t12+u22.*t22;
    u11=a11;u12=a12;u21=a21;u22=a22;
end
den=u22;
num1=zeros(ns,nmax);
for ii=1:ns
    u11=ones(1,nmax);u12=zeros(1,nmax);u21=u12;u22=u11;
    r=[s(1,:),0];
    Je=nnz(r>=re(ii));if Je==0, Je=1; end
    r(Je+1)=re(ii);
    for jj=1:Je
        [t11,t12,t21,t22]=multiUani(r(jj),r(jj+1),s(2:3,jj),nmax);
        a11=u11.*t11+u12.*t21;
        a12=u11.*t12+u12.*t22;
        a21=u21.*t11+u22.*t21;
        a22=u21.*t12+u22.*t22;
        u11=a11;u12=a12;u21=a21;u22=a22;
    end
    la=(repmat(s(1,Je)/re(ii),1,nmax).^((sqrt(1+4*(1:nmax).*((1:nmax)+1).*...
        repmat(s(3,Je)./s(2,Je),1,nmax))-1)/2))/re(ii)^2;
    for jj=(Je+1):size(s,2)
        nu=(sqrt(1+4*(1:nmax).*((1:nmax)+1).*repmat(s(3,jj-1)./s(2,jj-1),1,nmax))-1)/2;
        la=la.*(repmat(s(1,jj)./s(1,jj-1),1,nmax).^nu);
    end
    num1(ii,:)=u22.*la;
end
u=zeros(ns,nd);
for kk=1:nd
    p=di(kk,1:3);
    r0=norm(p);
    px=X*p'./re/r0;
    r=[s(1,:),0];J0=nnz(r>r0);r(J0)=r0;
    u11=ones(1,nmax);u12=zeros(1,nmax);u21=u12;u22=u11;
    for ii=J0:N
        [t11,t12,t21,t22]=multiUani(r(ii),r(ii+1),s(2:3,ii),nmax);
        a11=u11.*t11+u12.*t21;
        a12=u11.*t12+u12.*t22;
        a21=u21.*t11+u22.*t21;
        a22=u21.*t12+u22.*t22;
        u11=a11;u12=a12;u21=a21;u22=a22;
    end
    la=repmat(r0/s(1,J0),1,nmax).^((sqrt(1+4*(1:nmax).*((1:nmax)+1).*repmat(s(3,J0)./s(2,J0),1,nmax))-1)/2);
    for ii=(J0+1):N
        nu=(sqrt(1+4*(1:nmax).*((1:nmax)+1).*repmat(s(3,ii-1)./s(2,ii-1),1,nmax))-1)/2;
        la=la.*(repmat(s(1,ii-1)./s(1,ii),1,nmax).^nu);
    end
    la=(2*(1:nmax)+1).*la;
    Rn=la.*u12./den;
    dRn=la.*u22./den/s(2,J0);
    dRn2=s(3,J0)/s(2,J0)/r0^2*(1:nmax).*((1:nmax)+1).*Rn-2/r0*dRn; % nuevo
    P=zeros(ns,nmax+1);
    P(:,1)=ones(ns,1);
    P(:,2)=px;
    dP=zeros(ns,nmax);
    dP(:,1)=ones(ns,1);
    dP2=zeros(ns,nmax); % nuevo
    dP2(:,1)=3*ones(ns,1);
    for n=2:nmax
        P(:,n+1)=((2*n-1)*px.*P(:,n)-(n-1)*P(:,n-1))/n;
        dP(:,n)=n*P(:,n)+px.*dP(:,n-1);                  % Eq. 59
        dP2(:,n)=(n+2)*dP(:,n)+px.*dP2(:,n-1); % Correction from de Munck's paper (eq. 65): is (n+1)P'+P'' instead of nP'+xP''
    end
    P=P(:,2:nmax+1);
    dP2=[zeros(ns,1),dP2(:,1:nmax-1)];
    S0=sum(repmat(Rn,ns,1).*num1.*dP,2)/r0;          % Eq. 58
    S1=sum(repmat(dRn,ns,1).*num1.*P,2);             % Eq. 57 
 
    if size(di,2)==4 % monopole
        
        u(:,kk)=di(kk,4)*sum(repmat(Rn,ns,1).*num1.*P,2)/(4*pi);  
        
    elseif size(di,2)==6 % dipole
        
        q=di(kk,4:6);
        u(:,kk)=(q*p'/r0*(S1-px.*S0)+X*q'./re.*S0)/(4*pi);   % Eq. 56        
        
    elseif size(di,2)==9 % quadrupole
                
        S2=sum(repmat(dRn2,ns,1).*num1.*P,2);          
        S3=sum(repmat(dRn,ns,1).*num1.*dP,2)/r0;       
        S4=sum(repmat(Rn,ns,1).*num1.*dP2,2)/r0^2;     
        
        xo=p/r0;
        xe=X./repmat(re,1,3);
        
        Qd=eye(3);
        Qq=ivech(di(4:end),3);
        up=zeros(size(xe,1),3);u=up(:,1);
        for iq=1:3
            up=kron(xo,(Qd(iq,:)*xo')*(3*px.*S0/r0-S1/r0+S2-2*px.*S3+px.^2.*S4)+...
               dot(repmat(Qd(iq,:),size(xe,1),1),xe,2).*(-S0/r0+S3-px.*S4))+...
               kron(Qd(iq,:),(S1/r0-px.*S0/r0))+...        
               xe.*repmat((Qd(iq,:)*xo')*(-S0/r0+S3-px.*S4)+dot(repmat(Qd(iq,:),size(xe,1),1),xe,2).*S4,1,3);
            u=u+up*Qq(iq,:)';
        end
        u(:,kk)=sum(u,2)/(4*pi)/2;
        
    end
    
end


function [u11,u12,u21,u22]=multiUani(ra,rb,sig,nmax)

n=1:nmax;
nu=(sqrt(1+4*n.*(n+1)*sig(2)/sig(1))-1)/2;   % Eq. 14.
if rb~=0   
    c=1./(2*nu+1);     % Eq. 31 (with 29 y 30).
    d=(ra/rb).^(-2*(nu+1));
    u11=(nu.*d*ra/rb+nu+1).*c;
    u12=(rb-ra*d)/sig(1).*c;
    u21=sig(1)*nu.*(nu+1).*(rb-ra.*d)/(ra*rb).*c;
    u22=(rb/ra*nu+d.*(nu+1)).*c;
else
    c=1./(2*nu+1);     % Eq. 35.
    u11=0*c;
    u12=c/sig(1);
    u21=0*c;
    u22=nu/ra.*c;
end


function M=ivech(vect,num)

M=zeros(num);index=1;
for i=1:num
    for j=i:num
        M(i,j)=vect(index);
        M(j,i)=M(i,j);
        index=index+1;
    end
end
