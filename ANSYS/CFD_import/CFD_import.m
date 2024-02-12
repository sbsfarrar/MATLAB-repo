%% Copyright 2015 The MathWorks, Inc.

% This script is the main script for importing CFD simulation results
% and visualizing these. Furthermore it calculates
% momentum balance and drag force and evaluates the swirling strength.

%% import & format data
% import raw data
%be sure to have 'ANSYS' folder selected in present working directory
%(example filepath below)
path = pwd; 
% filepath = "C:\Users\StevenSummey\Documents\MATLAB\ANSYS\FFF-Setup-Output-Solution_ALL_ASCII";
filepath = [path '\FFF-Setup-Output-Solution_ALL_ASCII'];
CFD_raw = importdata(filepath); 

% replace spaces, then dashes with underscores in the variable names
variables_names_noSpaces = arrayfun(@(str) strtrim(str), CFD_raw.colheaders); 
variables_names = arrayfun(@(str) strrep(str, '-', '_'), variables_names_noSpaces); 
names_length = length(variables_names);
% build a table from the imported data
i = 1;
while i < names_length
    if string(variables_names{2}) == string(variables_names{i}) && i ~= 2
        indx1 = i;
        variables_names{i} = [variables_names{i} '1'];
    elseif string(variables_names{3}) == string(variables_names{i}) && i ~= 3
        indx2 = i;
        variables_names{i} = [variables_names{i} '1'];
    elseif string(variables_names{4}) == string(variables_names{i}) && i ~= 4
        indx3 = i;
        variables_names{i} = [variables_names{i} '1'];
    end
    i = i + 1;
end

% variables_names{indx1) = [];
% variables_names(indx2) = [];
% variables_names(indx3) = [];
    
CFD_data = array2table(CFD_raw.data, 'VariableNames', variables_names);

% get x, y, z min/max
x_min = min(CFD_data.x_coordinate);x_max = max(CFD_data.x_coordinate);
y_min = min(CFD_data.y_coordinate);y_max = max(CFD_data.y_coordinate);
z_min = min(CFD_data.z_coordinate);z_max = max(CFD_data.z_coordinate);

% get scattered interpolants
F_x = scatteredInterpolant(CFD_data.x_coordinate, CFD_data.y_coordinate, CFD_data.z_coordinate, CFD_data.x_velocity);
F_y = scatteredInterpolant(CFD_data.x_coordinate, CFD_data.y_coordinate, CFD_data.z_coordinate, CFD_data.y_velocity);
F_z = scatteredInterpolant(CFD_data.x_coordinate, CFD_data.y_coordinate, CFD_data.z_coordinate, CFD_data.z_velocity);
F_p = scatteredInterpolant(CFD_data.x_coordinate, CFD_data.y_coordinate, CFD_data.z_coordinate, CFD_data.total_pressure);

%% display a velocity slice quiver plot

% get meshgrid
spacing = 200;
y = linspace(y_min, y_max, spacing);
z = linspace(z_min, z_max, spacing);
[y_slice, z_slice] = meshgrid(y, z);

% interpolate slice data
slice_vx = reshape(F_x(zeros(spacing^2, 1), y_slice(:), z_slice(:)), [spacing spacing]);
slice_vy = reshape(F_y(zeros(spacing^2, 1), y_slice(:), z_slice(:)), [spacing spacing]);
slice_vz = reshape(F_z(zeros(spacing^2, 1), y_slice(:), z_slice(:)), [spacing spacing]);
slice_p  = reshape(F_p(zeros(spacing^2, 1), y_slice(:), z_slice(:)), [spacing spacing]);

% display slice
h_fig = figure(1);
set(h_fig, 'color', 'w');
hold on;
pcolor(y, z, sqrt(slice_vx.^2 + slice_vy.^2 + slice_vz.^2)); %set background
shading interp;
quiver_spacing = 4;
quiver_scale   = 3;

% draw quiver plot
h_quiv = quiver(y(1:quiver_spacing:end), z(1:quiver_spacing:end), ... 
                slice_vy(1:quiver_spacing:end, 1:quiver_spacing:end),...
                slice_vz(1:quiver_spacing:end, 1:quiver_spacing:end),...
                quiver_scale);

% set figure propertues
h_quiv.Color = [1 0 0];
set(gca, 'FontSize', 20);
xlabel('y [m]', 'FontSize', 30);
ylabel('z [m]', 'FontSize', 30);
cbar = colorbar;
cbar.Label.String = 'Velocity magnitude [m/s]';

%% effect momentum balance to compute drag and lift

% define inlet/outlet
y_outlet = -0.05;
y_inlet = 0.06;
x_lin = linspace(-0.06, 0.06, spacing);
z = linspace(-0.04, 0.04, spacing);
[x_slice, z_slice] = meshgrid(x_lin, z);

% inlet 
inlet_vx = F_x(x_slice(:), y_inlet*ones(spacing^2, 1), z_slice(:));
inlet_vy = F_y(x_slice(:), y_inlet*ones(spacing^2, 1), z_slice(:));
inlet_vz = F_z(x_slice(:), y_inlet*ones(spacing^2, 1), z_slice(:));
inlet_p  = F_p(x_slice(:), y_inlet*ones(spacing^2, 1), z_slice(:));

% outlet
outlet_vx = F_x(x_slice(:), y_outlet*ones(spacing^2, 1), z_slice(:));
outlet_vy = F_y(x_slice(:), y_outlet*ones(spacing^2, 1), z_slice(:));
outlet_vz = F_z(x_slice(:), y_outlet*ones(spacing^2, 1), z_slice(:));
outlet_p  = F_p(x_slice(:), y_outlet*ones(spacing^2, 1), z_slice(:));

surface = 0.12 * 0.08; % (m^2)
momentum_drag = mean(CFD_data.density)*(inlet_vy.*inlet_vy - outlet_vy.*outlet_vy);
pressure_drag = inlet_p - outlet_p;
drag = mean(momentum_drag + pressure_drag) * surface;
lift = mean(max(CFD_data.density)*(inlet_vy.*inlet_vz - outlet_vy.*outlet_vz)) * surface;

% add drag vector 
figure(h_fig);
arrow = annotation('arrow');
arrow.Color = 'green';
arrow.Position = [0.56 0.52 -drag/(drag + lift)/6 -lift/(drag + lift)/6];
arrow.LineWidth = 2.5;
arrow.HeadLength = 25;
arrow.HeadWidth = 15;

%% compute swirling strength
CFD_data.lambda_ci = arrayfun(@computeSwStr, CFD_data.dx_velocity_dx, CFD_data.dx_velocity_dy, CFD_data.dx_velocity_dz,...
                                             CFD_data.dy_velocity_dx, CFD_data.dy_velocity_dy, CFD_data.dy_velocity_dz,... 
                                             CFD_data.dz_velocity_dx, CFD_data.dz_velocity_dy, CFD_data.dz_velocity_dz);
                                
%% image plot of lambda Ci
% get scattered interpolant
F_lci = scatteredInterpolant(CFD_data.x_coordinate, CFD_data.y_coordinate, CFD_data.z_coordinate, CFD_data.lambda_ci);

spacing = 400;
x = linspace(x_min + 0.01, x_max - 0.01, spacing);
z = linspace(z_min + 0.01, z_max - 0.01, spacing);
y_pos = -0.085;
[x_slice, z_slice] = meshgrid(x, z);

% interpolate slice data
slice_lci = reshape(F_lci(x_slice(:), y_pos*ones(spacing^2, 1), z_slice(:)), [spacing spacing]);
figure(2);
pcolor(x, z, imfilter(slice_lci, fspecial('disk', 10)));

% set figure properties
set(gcf, 'Color', 'w');
shading interp;
cbar = colorbar;
cbar.Label.String = '\lambda_{ci} [s^{-1}]';
cbar.Label.FontSize = 20;
set(gca, 'FontSize', 20);
xlabel('y [m]', 'FontSize', 30);
ylabel('z [m]', 'FontSize', 30);




