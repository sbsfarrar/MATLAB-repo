% Example 2: solid model
clc; clear vars; close all;

% import the model
model = ANSYSimport('model2.txt');

% get mesh data for surface parts only
elements = model.surf.flist(:,1:4);
nodes    = model.surf.nlist(:,2:4);

% get results for surface parts only
Sx  = model.surf.stress(:,2);
Sy  = model.surf.stress(:,3);
Sz  = model.surf.stress(:,4);
Sxy = model.surf.stress(:,5);
Syz = model.surf.stress(:,6);
Sxz = model.surf.stress(:,7);
Ux  = model.surf.disp(:,2);
Uy  = model.surf.disp(:,3);
Uz  = model.surf.disp(:,4);

% open figure and plot
figure('name','Solid model showing displacements Ux')
patch('Faces',elements,'Vertices',nodes,'facevertexcdata',Ux,'facecolor','interp','edgecolor','k');
axis equal
axis vis3d
colormap(jet(9))
% axis off
xlim([-200 200])
ylim([-200 200])
zlim([-200 200])
view(45,45)
camlight

dragzoom