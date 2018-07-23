function [b,bv,bs]=femeg_indep_analyt(p,t,D,r0,q)


bs=-femeg_source_s(p,t,r0,q);
bv=-femeg_source_v(p,t,D,r0,q);
b=bv+bs;

end

function bs=femeg_source_s(p,t,r0,q)

% Find surface triangles
tri=surftri(p,t);

p1=p(tri(:,1),:);p2=p(tri(:,2),:);p3=p(tri(:,3),:);

% normal
n1=cross(p2-p1,p3-p1,2); % (53)
twoL=sqrt(sum(n1.^2,2)); % (54)
wn=n1./repmat(twoL,1,3); % (54)

l1=sqrt(sum((p3-p2).^2,2));
l2=sqrt(sum((p3-p1).^2,2));
l3=sqrt(sum((p1-p2).^2,2)); % (55)

un=(p2-p1)./repmat(l3,1,3);vn=cross(wn,un,2); % (56)

u3=dot(p3-p1,un,2);v3=twoL./l3; % (58)

u0=dot(un,[r0(1)-p1(:,1),r0(2)-p1(:,2),r0(3)-p1(:,3)],2);
v0=dot(vn,[r0(1)-p1(:,1),r0(2)-p1(:,2),r0(3)-p1(:,3)],2);
w0=-dot(wn,[r0(1)-p1(:,1),r0(2)-p1(:,2),r0(3)-p1(:,3)],2); % (59)


t0_1=(v0.*(u3-l3)+v3.*(l3-u0))./l1;
t0_2=(u0.*v3-v0.*u3)./l2;
t0_3=v0; % (4)

% (3)
sm_1=-((l3-u0).*(l3-u3)+v0.*v3)./l1;
sp_1=((u3-u0).*(u3-l3)+v3.*(v3-v0))./l1;
sm_2=-(u3.*(u3-u0)+v3.*(v3-v0))./l2;
sp_2=(u0.*u3+v0.*v3)./l2;
sm_3=-u0;
sp_3=l3-u0;

R0_1=sqrt(t0_1.^2+w0.^2);
R0_2=sqrt(t0_2.^2+w0.^2);
R0_3=sqrt(t0_3.^2+w0.^2);

Rp_1=sqrt(t0_1.^2+w0.^2+sp_1.^2);
Rp_2=sqrt(t0_2.^2+w0.^2+sp_2.^2);
Rp_3=sqrt(t0_3.^2+w0.^2+sp_3.^2);
Rm_1=sqrt(t0_1.^2+w0.^2+sm_1.^2);
Rm_2=sqrt(t0_2.^2+w0.^2+sm_2.^2);
Rm_3=sqrt(t0_3.^2+w0.^2+sm_3.^2);

bet_1=atan(t0_1.*sp_1./(R0_1.^2+abs(w0).*Rp_1))-atan(t0_1.*sm_1./(R0_1.^2+abs(w0).*Rm_1));
bet_2=atan(t0_2.*sp_2./(R0_2.^2+abs(w0).*Rp_2))-atan(t0_2.*sm_2./(R0_2.^2+abs(w0).*Rm_2));
bet_3=atan(t0_3.*sp_3./(R0_3.^2+abs(w0).*Rp_3))-atan(t0_3.*sm_3./(R0_3.^2+abs(w0).*Rm_3));
bet_t=bet_1+bet_2+bet_3;

I_w=bet_t./abs(w0);

sn_1=(p3-p2)./repmat(l1,1,3);sn_2=(p1-p3)./repmat(l2,1,3);sn_3=(p2-p1)./repmat(l3,1,3); % (55)
mn_1=cross(sn_1,wn,2);mn_2=cross(sn_2,wn,2);mn_3=cross(sn_3,wn,2); % (55)
f2_1=log((Rp_1+sp_1)./(Rm_1+sm_1));f2_2=log((Rp_2+sp_2)./(Rm_2+sm_2));f2_3=log((Rp_3+sp_3)./(Rm_3+sm_3)); % (12)

Rs_1=1./R0_1.^2.*(sp_1./Rp_1-sm_1./Rm_1);
Rs_2=1./R0_2.^2.*(sp_2./Rp_2-sm_2./Rm_2);
Rs_3=1./R0_3.^2.*(sp_3./Rp_3-sm_3./Rm_3);

% I0
Ip1=mn_1.*repmat(Rs_1,1,3)+mn_2.*repmat(Rs_2,1,3)+mn_3.*repmat(Rs_3,1,3);
Ip2=t0_1.*Rs_1+t0_2.*Rs_2+t0_3.*Rs_3;
I0=Ip1.*repmat(w0,1,3)-wn.*repmat(Ip2,1,3);
I0=q(1)*I0(:,1)+q(2)*I0(:,2)+q(3)*I0(:,3);

% I_u and I_v
st_1=sn_1.*repmat(1./Rm_1-1./Rp_1,1,3)+repmat(t0_1.*Rs_1,1,3).*mn_1;
st_2=sn_2.*repmat(1./Rm_2-1./Rp_2,1,3)+repmat(t0_2.*Rs_2,1,3).*mn_2;
st_3=sn_3.*repmat(1./Rm_3-1./Rp_3,1,3)+repmat(t0_3.*Rs_3,1,3).*mn_3;

mw=mn_1.*repmat(Rs_1.*w0.^2-f2_1,1,3)+mn_2.*repmat(Rs_2.*w0.^2-f2_2,1,3)+mn_3.*repmat(Rs_3.*w0.^2-f2_3,1,3);
I_u=repmat(w0,1,3).*(mn_1.*repmat(dot(un,st_1,2),1,3)+mn_2.*repmat(dot(un,st_2,2),1,3)+mn_3.*repmat(dot(un,st_3,2),1,3)-repmat(I_w,1,3).*un)+wn.*repmat(dot(un,mw,2),1,3);
I_u=q(1)*I_u(:,1)+q(2)*I_u(:,2)+q(3)*I_u(:,3);
I_v=repmat(w0,1,3).*(mn_1.*repmat(dot(vn,st_1,2),1,3)+mn_2.*repmat(dot(vn,st_2,2),1,3)+mn_3.*repmat(dot(vn,st_3,2),1,3)-repmat(I_w,1,3).*vn)+wn.*repmat(dot(vn,mw,2),1,3);
I_v=q(1)*I_v(:,1)+q(2)*I_v(:,2)+q(3)*I_v(:,3);

m12=-1./l3;m13=(u3./l3-1)./v3;m22=1./l3;m23=-u3./l3./v3;m33=1./v3;

I0_x=m12.*u0+m13.*v0+1;
I0_y=m22.*u0+m23.*v0;
I0_z=m33.*v0;

Iab_x=I_u.*m12+I_v.*m13;
Iab_y=I_u.*m22+I_v.*m23;
Iab_z=I_v.*m33;

I=[I0.*I0_x+Iab_x,I0.*I0_y+Iab_y,I0.*I0_z+Iab_z];


% Assemble vector
bs=zeros(size(p,1),1);
for kk=1:size(I,2)
    bs=bs+accumarray(tri(:,kk),I(:,kk)/4/pi,[size(p,1),1]);
end




end


function bv=femeg_source_v(p,t,D,r0,q)


Tpos=tsearchn(p,t(:,1:4),r0);
Dinf=D(Tpos(1),1);D(:,[1,4,6])=D(:,[1,4,6])-Dinf;


Ind_cn=sum(abs(D),2)>1e-12;

tri=[t(Ind_cn,[1,3,2]);t(Ind_cn,[1,2,4]);t(Ind_cn,[1,4,3]);t(Ind_cn,[3,4,2])];

[~,IA,IC]=unique(sort(tri,2),'rows');

Ind_sign=-ones(size(tri,1),1);Ind_sign(IA)=1;

Int_fa=int_rm1(p,tri(IA,:),r0,q);

I1=sum(reshape(Int_fa(IC,1).*Ind_sign,[],4),2);
I2=sum(reshape(Int_fa(IC,2).*Ind_sign,[],4),2);
I3=sum(reshape(Int_fa(IC,3).*Ind_sign,[],4),2);

% Compute volume coordintes
[b,c,d,V,~]=femeg_vol_coord(p,t(Ind_cn,1:4));

bev=[I1.*(b(:,1).*D(Ind_cn,1) + c(:,1).*D(Ind_cn,2) + d(:,1).*D(Ind_cn,3)) + I2.*(b(:,1).*D(Ind_cn,2) + c(:,1).*D(Ind_cn,4) + d(:,1).*D(Ind_cn,5)) + I3.*(b(:,1).*D(Ind_cn,3) + c(:,1).*D(Ind_cn,5) + d(:,1).*D(Ind_cn,6)),...
     I1.*(b(:,2).*D(Ind_cn,1) + c(:,2).*D(Ind_cn,2) + d(:,2).*D(Ind_cn,3)) + I2.*(b(:,2).*D(Ind_cn,2) + c(:,2).*D(Ind_cn,4) + d(:,2).*D(Ind_cn,5)) + I3.*(b(:,2).*D(Ind_cn,3) + c(:,2).*D(Ind_cn,5) + d(:,2).*D(Ind_cn,6)),...
     I1.*(b(:,3).*D(Ind_cn,1) + c(:,3).*D(Ind_cn,2) + d(:,3).*D(Ind_cn,3)) + I2.*(b(:,3).*D(Ind_cn,2) + c(:,3).*D(Ind_cn,4) + d(:,3).*D(Ind_cn,5)) + I3.*(b(:,3).*D(Ind_cn,3) + c(:,3).*D(Ind_cn,5) + d(:,3).*D(Ind_cn,6)),...
     I1.*(b(:,4).*D(Ind_cn,1) + c(:,4).*D(Ind_cn,2) + d(:,4).*D(Ind_cn,3)) + I2.*(b(:,4).*D(Ind_cn,2) + c(:,4).*D(Ind_cn,4) + d(:,4).*D(Ind_cn,5)) + I3.*(b(:,4).*D(Ind_cn,3) + c(:,4).*D(Ind_cn,5) + d(:,4).*D(Ind_cn,6))];
 

% Assemble vector
bv=zeros(size(p,1),1);
for kk=1:size(bev,2)
    bv=bv+accumarray(t(Ind_cn,kk),bev(:,kk)./V/6/4/pi/Dinf,[size(p,1),1]);
end

end



function I=int_rm1(p,t,r0,q)

p1=p(t(:,1),:);p2=p(t(:,2),:);p3=p(t(:,3),:);

% normal
n1=cross(p2-p1,p3-p1,2); % (53)
twoL=sqrt(sum(n1.^2,2)); % (54)
wn=n1./repmat(twoL,1,3); % (54)

l1=sqrt(sum((p3-p2).^2,2));
l2=sqrt(sum((p3-p1).^2,2));
l3=sqrt(sum((p1-p2).^2,2)); % (55)

un=(p2-p1)./repmat(l3,1,3);vn=cross(wn,un,2); % (56)

u3=dot(p3-p1,un,2);v3=twoL./l3; % (58)

u0=dot(un,[r0(1)-p1(:,1),r0(2)-p1(:,2),r0(3)-p1(:,3)],2);
v0=dot(vn,[r0(1)-p1(:,1),r0(2)-p1(:,2),r0(3)-p1(:,3)],2);
w0=-dot(wn,[r0(1)-p1(:,1),r0(2)-p1(:,2),r0(3)-p1(:,3)],2); % (59)

t0_1=(v0.*(u3-l3)+v3.*(l3-u0))./l1;
t0_2=(u0.*v3-v0.*u3)./l2;
t0_3=v0; % (4)

% (3)
sm_1=-((l3-u0).*(l3-u3)+v0.*v3)./l1;
sp_1=((u3-u0).*(u3-l3)+v3.*(v3-v0))./l1;
sm_2=-(u3.*(u3-u0)+v3.*(v3-v0))./l2;
sp_2=(u0.*u3+v0.*v3)./l2;
sm_3=-u0;
sp_3=l3-u0;

R0_1=sqrt(t0_1.^2+w0.^2);
R0_2=sqrt(t0_2.^2+w0.^2);
R0_3=sqrt(t0_3.^2+w0.^2);

Rp_1=sqrt(t0_1.^2+w0.^2+sp_1.^2);
Rp_2=sqrt(t0_2.^2+w0.^2+sp_2.^2);
Rp_3=sqrt(t0_3.^2+w0.^2+sp_3.^2);
Rm_1=sqrt(t0_1.^2+w0.^2+sm_1.^2);
Rm_2=sqrt(t0_2.^2+w0.^2+sm_2.^2);
Rm_3=sqrt(t0_3.^2+w0.^2+sm_3.^2);

bet_1=atan(t0_1.*sp_1./(R0_1.^2+abs(w0).*Rp_1))-atan(t0_1.*sm_1./(R0_1.^2+abs(w0).*Rm_1));
bet_2=atan(t0_2.*sp_2./(R0_2.^2+abs(w0).*Rp_2))-atan(t0_2.*sm_2./(R0_2.^2+abs(w0).*Rm_2));
bet_3=atan(t0_3.*sp_3./(R0_3.^2+abs(w0).*Rp_3))-atan(t0_3.*sm_3./(R0_3.^2+abs(w0).*Rm_3));

bet_t=bet_1+bet_2+bet_3;

I_w=bet_t./abs(w0);

%
sn_1=(p3-p2)./repmat(l1,1,3);
sn_2=(p1-p3)./repmat(l2,1,3);
sn_3=(p2-p1)./repmat(l3,1,3); % (55)

mn_1=cross(sn_1,wn,2);
mn_2=cross(sn_2,wn,2);
mn_3=cross(sn_3,wn,2); % (55)

f2_1=log((Rp_1+sp_1)./(Rm_1+sm_1));
f2_2=log((Rp_2+sp_2)./(Rm_2+sm_2));
f2_3=log((Rp_3+sp_3)./(Rm_3+sm_3)); % (12)

I_s=mn_1.*repmat(f2_1,1,3)+mn_2.*repmat(f2_2,1,3)+mn_3.*repmat(f2_3,1,3);

I=I_s-repmat(I_w.*w0,1,3).*wn;

%
I=-wn.*repmat(I(:,1)*q(1)+I(:,2)*q(2)+I(:,3)*q(3),1,3);

end





