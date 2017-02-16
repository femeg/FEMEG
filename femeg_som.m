function [pso,tso,inda] = femeg_som(p,t)
% FEMEG_SOM computes the second order mesh from the first order mesh  
%
% [pso,tso,inda] = FEMEG_SOM(p,t) computes the second order mesh with nodes
% pso and elements tso from the first order mesh with nodes p and elements
% t. The relation between these meshes is given by inda: pso(inda,:)=p;
% Also, p(t(:,ii),jj)=pso(tso(:,ii),jj), for all ii=1,2,3,4 and jj=1,2,3.
%
%
% Inputs:
%    p: nodes defining the first ordertetrahedral mesh (size: Np x 3).
%    t: elements defining the first order tetrahedral mesh (size: Nt x 4(+1)).
%
% Outputs:
%    pso: nodes defining the second order tetrahedral mesh (size: Npso x 3).
%    tso: elements defining the second order tetrahedral mesh (size: Ntso x 10(+1)).
%    inda: vector for relating pso and p (size: Np x 1).

%   This file is part of the FEMEG toolbox.
%   Author: Leandro Beltrachini <BeltrachiniL at cardiff.ac.uk>


p1=p(t(:,1),:);p5=p(t(:,2),:);p8=p(t(:,3),:);p10=p(t(:,4),:);
p2=(p1+p5)/2;p6=(p8+p5)/2;p9=(p10+p8)/2;p4=(p1+p10)/2;p7=(p10+p5)/2;p3=(p1+p8)/2;
psv=[p1;p2;p3;p4;p5;p6;p7;p8;p9;p10];
clear p2 p3 p4 p5 p6 p7 p8 p9 p10

dt=size(p1,1);
t1=1:dt;t2=dt+1:2*dt;t3=2*dt+1:3*dt;t4=3*dt+1:4*dt;t5=4*dt+1:5*dt;t6=5*dt+1:6*dt;
t7=6*dt+1:7*dt;t8=7*dt+1:8*dt;t9=8*dt+1:9*dt;t10=9*dt+1:10*dt;

[pso,~,IC]=unique(psv,'rows');clear psv
tso=[IC(t1),IC(t2),IC(t3),IC(t4),IC(t5),IC(t6),IC(t7),IC(t8),IC(t9),IC(t10)];

if size(t,2)==5,tso=[tso,t(:,end)];end

tso=[tso(:,1),tso(:,5),tso(:,8),tso(:,10),tso(:,2),tso(:,3),tso(:,4),tso(:,6),tso(:,7),tso(:,9),tso(:,11)];


if nargout==3
    hh=sqrt(min(sum((p(t(:,1),:)-p(t(:,2),:)).^2,2)));
    in1=10^(round(log10(hh))-2);

    [~,indb]=ismember(round(pso/in1),round(p/in1),'rows');
    indd=find(indb);
    indb(indb==0)=[];
    [~,I]=sort(indb);

    inda=indd(I);
end


