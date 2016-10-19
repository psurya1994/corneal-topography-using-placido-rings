function [r_rings, no_of_rings] = psa(r_cell, center, r_contour0, r_max, bool_plot)

%   r_contour0 -> shape of innermost ring

    % Initializing the parameters
    state = 1; 
    r_contour = r_contour0;
    no_of_rings = 1; % keep increasing as it goes inside the function
    r_rings{no_of_rings} = r_contour0;
    
    % Constructing r_matrix
    r_matrix = zeros(length(r_contour0), r_max);
    for i = 1:length(r_contour0)
        r_matrix(i, r_cell{i}) = 1;
    end

    r_contour = r_contour + 1;
    
    % Loop increasing value of a and b; end when max values of either of them
    % is crossed
    if(bool_plot == 1)
        figure
        for i = 1:length(r_cell)
            if(i == 1) 
                hold off
            end
            plot(i, r_cell{i},'.')
            hold on
        end
        xlabel('theta_k'), ylabel('r (distance from center)')
        hold on
    end
    
    while(max(r_contour) <= r_max)
        
        if(bool_plot == 1)
            plot(1:length(r_cell), r_contour, 'r-');
            pause(0.001);
        end
        
        disp(['state = ' num2str(state)])
        if(mod(state, 3) == 1) % growing to detect white
            if (~is_contour_on_dot(r_contour, r_matrix))
                r_contour = r_contour + 1;
            else
                r_start = r_contour;
                state = state + 1;
            end

        elseif (mod(state, 3) == 2) % detected white, finding black
            % See if all pixels on ellipse are 0
            % if (all not 0) -> radius increase
            % else -> stop
            if (is_contour_on_dot(r_contour, r_matrix))
                r_contour = r_contour + 1;
            else
                
                state = state - 1;
                r_end = r_contour;

                % extracting points between r_contour_start and r_contour_end
                % and storing in ring_x
                ring_x = size(r_contour);
                for i = 1:size(r_contour, 2)
                    if(isnan(r_start(i)) || isnan(r_end(i)))
                        ring_x(i) = nan;
                        continue
                    end
                    ind = find(r_matrix(i, r_start(i):r_end(i)) == 1, 1);
                    if(~isempty(ind))
                        ring_x(i) = r_start(i) + ...
                                find(r_matrix(i, r_start(i):r_end(i)) == 1);
                    else
                        ring_x(i) = nan;
                    end
                end

                % fit all internal points into a new contour
                r_contour = ring_x;
                no_of_rings = no_of_rings + 1;
                r_rings{no_of_rings} = ring_x - 1; % -1: bug compensation, needs fixing later

                r_contour = r_contour + 2; % 2 to avoid erroneous points situtated just next
            end 
        end
    end
end
