function [P, T] = femeg_sphere( r, N, m)
%FEMEG_SPHERE provides a spherical surface mesh where to test algorithms. 
%
%   [P, T] = femeg_sphere ( r, N) provides a set of planar triangles on a
%   spherical surface of radius r. The number of triangles is related with
%   N and m in the following way: Nt=(20*4^N)*(m+1) (m=0 by default).
%
%   [P, T] = femeg_sphere ( r, N, m) specifies m (m=0 or m=1).

%   This file is part of the FEMEG toolbox.
%   Author: Nicolas von Ellenrieder and Leandro Beltrachini <BeltrachiniL at cardiff.ac.uk>

if nargin==2 || m==0
    
    v=[1 -1 -1];
    tau=max(roots(v));
    P=[0  tau  1; 0  tau  -1; 0  -tau 1; 0  -tau -1; 1  0 tau; -1 0 tau;
            1  0 -tau; -1 0 -tau; tau 1 0; tau -1 0; -tau 1 0; -tau -1 0]/2;
    T=[2,1,9;9,7,2;11,1,2;2,8,11;6,1,11;11,12,6;
           5,1,6;6,3,5;9,1,5;5,10,9;10,4,7;7,9,10;
           7,4,8;8,2,7;8,4,12;12,11,8;12,4,3;3,6,12;3,4,10;10,5,3];
    clear tau v

else
    o=48*pi/180;
    P=[0,0,1;cos(o)*[1,0;cos(2*pi/5),sin(2*pi/5);cos(4*pi/5),sin(4*pi/5);
        cos(6*pi/5),sin(6*pi/5);cos(8*pi/5),sin(8*pi/5)],sin(o)*ones(5,1);
        1,0,0;cos(pi/5),sin(pi/5),0;cos(2*pi/5),sin(2*pi/5),0;
        cos(3*pi/5),sin(3*pi/5),0;cos(4*pi/5),sin(4*pi/5),0;-1,0,0;
        cos(6*pi/5),sin(6*pi/5),0;cos(7*pi/5),sin(7*pi/5),0;
        cos(8*pi/5),sin(8*pi/5),0;cos(9*pi/5),sin(9*pi/5),0;
        cos(o)*[1,0;cos(2*pi/5),sin(2*pi/5);cos(4*pi/5),sin(4*pi/5);
        cos(6*pi/5),sin(6*pi/5);cos(8*pi/5),sin(8*pi/5)],-sin(o)*ones(5,1);0,0,-1];
    T=[1,2,3;1,3,4;1,4,5;1,5,6;1,6,2;2,7,8;2,8,3;3,8,9;3,9,10;4,3,10;4,10,11;
        4,11,12;5,4,12;5,12,13;5,13,14;6,5,14;6,14,15;6,15,16;2,6,16;2,16,7;
        8,7,17;8,17,18;8,18,9;10,9,18;10,18,19;10,19,11;12,11,19;12,19,20;12,20,13;
        14,13,20;14,20,21;14,21,15;16,15,21;16,21,17;16,17,7;22,18,17;22,19,18;
        22,20,19;22,21,20;22,17,21];
    clear o
end

for j=1:N
	[P,T]=triang4(P,T);		
	P=P.*repmat(r./sqrt(sum(P.^2,2)),1,3);	
end



function [P,T]=triang4(P,T)

% This functions multiplies by 4 the number of triangles on a surface, adding a node in the midpoint of every triangle side

m=size(T,1);
n=size(P,1);
A=T;A(:,:,2)=[T(:,2:3),T(:,1)];A=sort(A,3);A=A(:,:,2)+n*A(:,:,1);
a=zeros(n^2,1);a(A(:))=1;a=find(a);
for ii=1:m
    p1=(P(T(ii,1),:)+P(T(ii,2),:))/2;a1=n+find(a==A(ii,1));P(a1,:)=p1;
    p2=(P(T(ii,2),:)+P(T(ii,3),:))/2;a2=n+find(a==A(ii,2));P(a2,:)=p2;
    p3=(P(T(ii,3),:)+P(T(ii,1),:))/2;a3=n+find(a==A(ii,3));P(a3,:)=p3;
    Tn(4*ii,:)=[a1,T(ii,2),a2];
    Tn(4*ii+1,:)=[a2,a3,a1];
    Tn(4*ii+2,:)=[a3,a2,T(ii,3)];
    Tn(4*ii+3,:)=[T(ii,1),a1,a3];
end
T=Tn(4:end,:);
