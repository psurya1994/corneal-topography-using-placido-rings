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
% 
% img = imread('data/test-1.png'); 
% center = [385, 427]; a_max = 300; b_max = 300;
% a0 = 5; b0 = 6;

% test - image 3
img = imread('data/test-4.png'); center = [385, 427];

figure, imshow(img);
%% 
segments = 36;
[points, r_cell] = intersection_points(img, center, segments, 0);

r_contour0 = zeros(1, segments);
for i = 1:segments
    r_contour0(i) = min(r_cell{i});
end
r_max = 180;

[r_rings, no_of_rings] = psa(r_cell, center, r_contour0, r_max, 0);

% doesn't work for some values of segments - code needs to be improved
% change segments value if code breaks
% rings_size = size(points{1}, 1);
% for i = 2:size(points, 1)
%     if(rings_size ~= size(points{i}, 1))
%         disp('ERROR. Please choose different value for segments.')
%         return
%     end
% end

%% Finding the radius of curvature
% make sure the calibrated surface paramters are stored as sf is in your 
% workspace

% r = points2r(points, center); % r (segment, ring_number)

for i = 1:no_of_rings
    r(:, i) = r_rings{i};
end

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
% figure, polarplot3d(rho_c');
figure, polarplot3d(rho_c(end:-1:1,:)');
caxis([7.5 13])
view([0 90]);
title('Radius of curvature plot')
colorbar

%% Plotting the diopter value
figure, polarplot3d(337./rho_c(end:-1:1,:)');
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