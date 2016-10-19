function [img_lb, ring_no] = esa(img, center, a0, b0, a_max, b_max, bool_plot) 
%ESA  an algorithm to detect and label rings in a placido ring
%image.
% 
% The algorithm has been explained in Espinosa's paper:
% https://rua.ua.es/dspace/bitstream/10045/25887/1/espinosa_et_al.pdf
% 
% Ellipse fit has been used from:
% Ohad Gal's fit_ellipse on MATLAB's file exchange:
% https://in.mathworks.com/matlabcentral/fileexchange/3215-fit-ellipse
% 
% The typical call of the function looks like this:
% [img_lb, ring_no] = ESA(img, center, a0, b0, a_max, b_max, bool_plot)
% 
% Input:
% img -> input binary image with the placido rings in white
% center -> center of the placido disk image in coordinates [x, y]
% a0 -> starting size of vertical axis of ellipse (in pixels)
% b0 -> starting size of horizontal axis of ellipse (in pixels)
% a_max -> maximum values of ellipse vertical axis (in pixels)
% b_max -> maximum values of ellipse horizontal axis (in pixels)
% boot_plot: '0' - no plots, '1' - plots each ellipse fit, '2' - plots the
% growing ellipse and ellipse fits
% 
% Output:
% img_lb -> output image with rings labelled with different numbers.
% imagesc can be used to display this image
% ring_no -> returns the number of rings in the image
% 
% Note: Tilt of the ellipse has not been used for updating ellipse
% parameters due to which the code might be ineffective at times.

% Initializing the parameters
a = a0;
b = b0;

% 'state' represents state of algorithm as it scans
% 
% state = 1 => ellipse scanning and increasing to find a white pixel
% 
% state = 2 => white pixel found, increasing in size until all elements on
% the perimeter of ellipse are black
% 
% state = 3 => all elements in perimeter are black; so fit an ellipse and
% update parameters accordingly

state = 1; 


img_lb = zeros(size(img));
ring_no = 0;

% Loop increasing value of a and b; end when max values of either of them
% is crossed
while(a < a_max || b < b_max)
    
    % Increment step size
    step_a = a/(a+b);
    step_b = b/(a+b);
    
    % Finding the points for the ellipse
    peri = 2 * pi * sqrt((a^2+b^2)/2);
    
    % any number more than 0.05 as multiplication factor works fine
    % this can be reduced to speed up process
    segments = round(peri * 0.1); 
    theta_step = (2*pi)/(segments);
    
    % Finding points using parametric method
    theta = 0 : theta_step : 2*pi-theta_step;
    x = center(1) + a * cos(theta);
    y = center(2) + b * sin(theta);

    % Converting the points into whole numbers (pixel values)
    % these pixels have redundant values that can be removed later
    x_px = round(x);
    y_px = round(y);

    % ploting the ellipse increasing in size
    if (bool_plot == 2)
        figure(1), imshow(img)
        hold on, plot(y_px, x_px,'r.')
    end
    
    if(mod(state, 3) == 1) % growing to detect white
        % See if there is a 1 detected on the ellipse
        % if (1 not detected) -> radius increase
        % else -> stop
        if (~is_pixel_on_image(img, x_px, y_px, 1))
            a = a + step_a;
            b = b + step_b;
        else
            a_start = a;
            b_start = b;
            state = state + 1;
        end    
        
    elseif (mod(state, 3) == 2) % detected white, finding black
        % See if all pixels on ellipse are 0
        % if (all not 0) -> radius increase
        % else -> stop
        if (is_pixel_on_image(img, x_px, y_px, 1))
            a = a + 2*a/(a+b);
            b = b + 2*b/(a+b);
        else
            a_end = a;
            b_end = b;
            
            % extract all internal points
            outer_ellipse = ellipse_mask(size(img), center, a_end, b_end, 0);
            inner_ellipse = ellipse_mask(size(img), center, a_start, b_start, 0);
            roi = (outer_ellipse & ~inner_ellipse) & (img);
            
            % fit all internal points into an ellipse
            stats = regionprops(roi, 'PixelList');
            if(length(stats)~=1)
                disp('Error: Region of interest- incorrectly calculated')
                figure, imshow(roi)
                beep
                return
            end
            
            if(bool_plot == 1)
                h = figure(1);
                imshow(img), hold on
            end
            
            try
                if(bool_plot == 1)
                    ellipse_t = fit_ellipse( stats.PixelList(:,1), stats.PixelList(:,2), h);
                else
                    ellipse_t = fit_ellipse( stats.PixelList(:,1), stats.PixelList(:,2));
                end
                ring_no = ring_no + 1;
                img_lb = img_lb + ring_no * roi;
            catch
                disp('Ellipse fit failed!')
                return
            end
            
            % Update parameters based on the ellipse found
            a = ellipse_t.b;
            b = ellipse_t.a;
            center(2) = ellipse_t.X0_in;
            center(1) = ellipse_t.Y0_in;
            state = state + 1;
            pause(1);
        end
    
    elseif(mod(state, 3)==0) % increasing new ellipse till all perim pixels black
        if (is_pixel_on_image(img, x_px, y_px, 1))
            a = a + step_a;
            b = b + step_b;
        else
            a_start = a;
            b_start = b;
            state = state + 1;
        end
    end
    
end

