function [V,h]=femeg_vis3d(p,t,d,expr)
% FEMEG_VIS3D plots scalar fields in 3d tetrahedral meshes
%
% [V,h] = FEMEG_VIS3D(p,t,d) plots the 3D first order tetrahedral mesh and
% colours it according to the scalar field d.
%
% [V,h] = FEMEG_VIS3D(p,t,d,expr) does the same but cutting the mesh
% according to expr.
%
% Inputs:
%    p: nodes defining the first ordertetrahedral mesh (size: Np x 3).
%    t: elements defining the first order tetrahedral mesh (size: Nt x 4(+1)).
%    d: scalar field map for colouring the mesh (size: Np x 1).
%    expr: expression for cutting the volume (string). Use "p(:,1)" for
%          referring to the "x-axis", "p(:,2)" for referring to the "y-axis",
%          and "p(:,3)" for referring to the "z-axis". E.g. expr='p(:,2)>0'.

%   This file is part of the FEMEG toolbox.
%   Author: Leandro Beltrachini <BeltrachiniL at cardiff.ac.uk>


tri1=surftri(p,t);
icol=.1*ones(1,3);
t=t(:,1:4);tor=t;

if size(t)<1,return,end


if nargin>3 && ~isempty(expr)
    incl=find(eval(expr));
    t=t(any(ismember(t,incl),2),:);
    tri1=tri1(any(ismember(tri1,incl),2),:);
    tri2=surftri(p,t);
    tri2=setdiff(tri2,tri1,'rows');
    h=trimesh(tri2,p(:,1),p(:,2),p(:,3),d(:));
    set(h,'facecolor',icol,'edgecolor','k');
    hold on
end
  
if size(d,1)==size(p,1)
    s.FaceVertexCData=d(:);s.Faces=tri1;
else
    ty=tsearchn(p,double(tor(:,1:4)),(p(tri2(:,1),:)+p(tri2(:,2),:)+p(tri2(:,3),:))/3);s.FaceVertexCData=d(ty);s.Faces=tri2;
end

s.Vertices=p;

V=patch(s);

if size(d(:),1)==size(p,1), shading interp; else shading flat; end

axis equal



function tri=surftri(p,t)
%SURFTRI Find surface triangles from tetrahedra mesh
%   TRI=SURFTRI(P,T)

%   Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.

% Form all faces, non-duplicates are surface triangles
faces=[t(:,[1,2,3]);
       t(:,[1,2,4]);
       t(:,[1,3,4]);
       t(:,[2,3,4])];
node4=[t(:,4);t(:,3);t(:,2);t(:,1)];
faces=sort(faces,2);
[~,ix,jx]=unique(faces,'rows');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
tri=faces(ix(qx),:);
node4=node4(ix(qx));

% Orientation
v1=p(tri(:,2),:)-p(tri(:,1),:);
v2=p(tri(:,3),:)-p(tri(:,1),:);
v3=p(node4,:)-p(tri(:,1),:);
ix=find(dot(cross(v1,v2,2),v3,2)>0);
tri(ix,[2,3])=tri(ix,[3,2]);

