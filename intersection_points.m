function [ points ] = intersection_points(img_lb, center, ring_no, segments, bool_plot)
%INTERSECTION_POINTS detects intersection points
%   Function returns points cell which is a {n, 1} cell where n represents
%   the number of segments. Inside each entry of the cell, you have the x
%   and y coordinates of all points lying on a specific segment.

    points = cell(segments,1); % filled with nan if point missing

    for ind = 1:ring_no % working with subimage with only ith ring
        sub_img = (img_lb == ind);
        i = 1;
        for theta = 0:360/segments:359.9
            theta_rad = theta*pi/180;
            search_point = center;
            r = 1; % step size for search point
            prev_spvalue = false; % previous search point value
            found = false;
            while(inside_bounds_sp(search_point, size(img_lb))==true)
                if(sub_img(search_point(1), search_point(2)) == ~prev_spvalue)
                    points{i} = [points{i}; search_point];
                    prev_spvalue = ~prev_spvalue;
                    found = true;
                    r = r + 1; % don't detect one pixel wide
                end
                % update search point
                search_point = center + ...
                            [floor(r*sin(theta_rad)), floor(r*cos(theta_rad))];
                r = r + 1;
            end
            if(found == false) % no color changes found in sub image
                points{i} = [points{i}; [nan nan]];
                points{i} = [points{i}; [nan nan]];
            end
            i = i + 1;
        end
    end

    % Plotting the dots
    if(bool_plot == true)
        figure, imshow(img_lb)
        for j = 1:length(points)
            hold on, plot(points{j}(:,2),points{j}(:,1), 'b.')
        end
    end
    
end

function in = inside_bounds_sp(point, img_size)
% INSIDE_BOUNDS_SP function returns if a search given point is inside the 
% bounds of the image.

    if(point(1)>0 && point(1)<img_size(1) ...
       && point(2)>0 && point(2)<img_size(2))
        in = true;
    else
        in = false;
    end
    
end