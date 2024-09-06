% Example 1: planar model
clc; clear vars; close all;

% import the model
model = ANSYSimport('model1.txt');

% get mesh data for plotting
elements      = model.raw.flist(:,1:4);
node_numbers  = model.raw.nlist(:,1);
nodes_x       = model.raw.nlist(:,2);
nodes_y       = model.raw.nlist(:,3);

% get results:
Sx  = model.raw.stress(:,2);
Sy  = model.raw.stress(:,3);
Sxy = model.raw.stress(:,5);
Ux  = model.raw.disp(:,2);
Uy  = model.raw.disp(:,3);


% vertices array
verts = [nodes_x nodes_y];

% same, including scaled displacements
d_scale = 5000;
verts_disp = [nodes_x+d_scale*Ux nodes_y+d_scale*Uy];

% make cell array of node number strings
for i = 1:length(node_numbers)
    n_str{i} = num2str(node_numbers(i));
end


% open figure and plot
figure('name','Planar model')

subplot(1,2,1) % model with node numbers
patch('Faces',elements,'Vertices',verts,'facecolor','c','edgecolor','k');
axis equal
text(nodes_x,nodes_y,n_str,'color','r','fontsize',8,'VerticalAlignment','bottom'); 


subplot(1,2,2) % stress in x-direction and deformed mesh
patch('Faces',elements,'Vertices',verts,'facevertexcdata',Sx,'facecolor','interp','edgecolor','none');
axis equal
colormap jet
patch('Faces',elements,'Vertices',verts_disp,'facevertexcdata',Sx,'facecolor','none','edgecolor','k');