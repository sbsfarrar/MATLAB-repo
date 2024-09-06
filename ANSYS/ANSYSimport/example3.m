% Example 3: shell model
clc; clear vars; close all;

% import the model
model = ANSYSimport('model3.txt');

% get mesh data for plotting
elements = model.raw.flist(:,1:4);
nodes    = model.raw.nlist(:,2:4);

% get results:
Sx  = model.raw.stress(:,2);
Sy  = model.raw.stress(:,3);
Sz  = model.raw.stress(:,4);
Sxy = model.raw.stress(:,5);
Syz = model.raw.stress(:,6);
Sxz = model.raw.stress(:,7);
Ux  = model.raw.disp(:,2);
Uy  = model.raw.disp(:,3);
Uz  = model.raw.disp(:,4);

% use displaced nodal coordinates
d_scale = 100;
nodes = nodes + d_scale*[Ux Uy Uz];

% open figure and plot
figure('name','Shell model, deformed, showing Sx')
patch('Faces',elements,'Vertices',nodes,'facevertexcdata',Sx,'facecolor','interp','edgecolor','k');
axis equal
axis vis3d
colormap(jet(11))
% axis off
view(210,30)
camlight

dragzoom