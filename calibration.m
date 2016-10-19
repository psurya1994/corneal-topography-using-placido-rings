%% Calibration using standard spheres
% 
% This fits a surface through the data points of your calibration spheres 
% and stores it in a variable `sf`. This step might take more than a minute
% to run.

tic
% Setting things up
clear all; clc; close all;


% Initialzing required parameters
spheres = 12;
rc = 8:0.5:13.5; % because the simulated rings have this radius
center = [385, 427]; 
a_max = 300; b_max = 300; 
a0 = 20; b0 = 20;
segments = 36;
r_spheres = cell(spheres, 1);
L = 18; % number of rings

r_contour0 = zeros(1, segments);
r_max = 300;

h = waitbar(0, 'Calibrating spheres...');

% Looping through each of the calibration sphere images
for i = 1:spheres
    % load image
    img = imread(['data/calib-' num2str(i) '.png']);
    
    % elliptical scanning
%     [img_lb, ring_no] = esa(img, center, a0, b0, a_max, b_max, 0);
%     waitbar((i-0.5) / spheres, h)
    
    % finding intersection points
%     points = intersection_points(img_lb, center, ring_no, segments);
    
    [points, r_cell] = intersection_points(img, center, segments, 0);
    r_contour0 = zeros(1, segments);
    r_max = 280;
    for m = 1:segments
        r_contour0(m) = min(r_cell{m});
    end
    [r_rings, no_of_rings] = psa(r_cell, center, r_contour0, r_max, 0);
    
    for m = 1:no_of_rings
        r(:, m) = r_rings{m};
    end

    % converting to right format
    [ r_spheres{i} ] = r;
    
    waitbar(i / spheres, h)
end
close(h)

%% Curve fitting to find the parameters of the fitting line


surface_l = [];
surface_rho = [];
surface_r = [];

figure
for l = 1:L
    r_data = [];
    rho_data = [];
    for j = 1:spheres
        r_data = [r_data r_spheres{j,1}(:,l)'];
        rho_data = [rho_data repmat(rc(j), 1, segments)];
    end
    p = polyfit(r_data, rho_data, 1);
    surface_l = [surface_l l*ones(1,length(r_data))]; %#ok<*AGROW>
    surface_rho = [surface_rho rho_data];
    surface_r = [surface_r r_data];
    hold on, plot(l, p(1), 'rx')
    hold on, plot(l, p(2), 'r.')
end
title('The fitting paramters')
xlabel('ring number'), ylabel('m_l or n_l')


%% Plotting calibration surface

% points plot
figure, plot3(surface_l, surface_rho, surface_r, 'ro')
title('Calibration surface plot')
xlabel('ring number (l)'), ylabel('rho (mm)'), zlabel('R (pixels)')

% surface fitting and plot
sf = fit([surface_l', surface_rho'],surface_r','poly21');
figure, plot(sf,[surface_l',surface_rho'],surface_r')
title('Calibration surface plot')
xlabel('ring number (l)'), ylabel('rho (mm)'), zlabel('R (pixels)')

disp('Calibration successful! The parameters of linear surface fit are:')
disp(sf)
toc