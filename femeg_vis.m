function V=femeg_vis(p,t,d)

s.Vertices=p;
s.Faces=t;
s.FaceVertexCData=d(:);
V=patch(s);
if size(d(:),1)==size(p,1), shading interp; end
if size(d(:),1)==size(t,1), shading flat; end

axis equal

