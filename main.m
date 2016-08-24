%% Generating Polar Corneal Topography 

close all
%% Loading the image whose radius of curvature must be found 
% Uncomment the image you'd like to load and comment others

% test image - 1: this is one of the input images on which the surface it
% fit; in fact this is a first sphere of 8mm radius of curvature. Hence,
% the result should have an approximate 8mm radius of curvature all over.
% 
% Comment/uncomment below lines for test image - 1.

% img = imread('data/calib-1.png'); 
% center = [385, 427]; a_max = 300; b_max = 300;
% a0 = 10; b0 = 10;

% test image - 2: this is a distorted input image that has been expanded in
% the horizontal direction. Hence, the result must give a higher radius of
% curvature in horizontal direction and less in vertical.
% 
% Comment/uncomment below lines for test image - 2.

img = imread('data/test-1.png'); 
center = [385, 427]; a_max = 300; b_max = 300;
a0 = 5; b0 = 6;

%% Elliptical scanning alogirthm for labelling image
[img_lb, ring_no] = esa(img, center, a0, b0, a_max, b_max, 1);
figure, imagesc(img_lb)

%% Detect intersection points with lines drawn from center
segments = 24; % should be same as when it's calibrated 
points = intersection_points(img_lb, center, ring_no, segments, 1);

% doesn't work for some values - code needs to be improved
%  if a non-working value of segements is selected, code breaks
rings_size = size(points{1}, 1);
for i = 2:size(points, 1)
    if(rings_size ~= size(points{i}, 1))
        disp('ERROR. Please choose different value for segments.')
        return
    end
end

%% Finding the radius of curvature
% make sure the calibrated surface paramters are stored as sf is in your 
% workspace

r = points2r(points, center); % r (segment, ring_number)
rho_c = size(r); % rho_c (segment, ring_number)
for i = 1:size(r, 1)
    for j = 1:size(r, 2)
        % first degree approx
        % rho_c(i, j) = (1/sf.p01) * r(i,j) - (sf.p10/sf.p01) * j - (sf.p00/sf.p01);
        % second degree in l, first degree in pho approx
        rho_c(i, j) = (r(i,j)-sf.p00-sf.p10*j-sf.p20*j^2)/(sf.p01+sf.p11*j);
    end
end

%% Plotting the curvature
figure, polarplot3d(rho_c');
caxis([7.5 13])
view([0 90]);
title('Radius of curvature plot')
colorbar

%% Plotting the diopter value
figure, polarplot3d(337./rho_c');
caxis([26 45])
view([0 90]);
title('Diopter plot')
colorbar

%% Plotting R's of all the rings
colors = ['y', 'm', 'c', 'r', 'g', 'b', 'k'];

figure
for i = 1:size(r,2)
    hold on, plot(r(:,i), colors(mod(i,length(colors))+1));
    legendInfo{i} = ['ring ' num2str(i)]; %#ok<SAGROW>
end
legend(legendInfo, 'Location', 'NorthEastOutside')
title('Distance of points on various rings from center')
xlabel('theta_k'), ylabel('R (pixels)')

%% Plotting rho's of all the rings

figure
for i = 1:size(r,2)
    hold on, plot(rho_c(:,i), colors(mod(i,length(colors))+1));
    legendInfo{i} = ['ring = ' num2str(i)]; %#ok<SAGROW>
end
legend(legendInfo, 'Location', 'NorthEastOutside')
title('Radius of curvature of points on various rings from center')
xlabel('theta_k'), ylabel('rho (mm)')