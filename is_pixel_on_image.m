function res = is_pixel_on_image(img, x_px, y_px, search_val)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    res = false;
    for i = 1:length(x_px)
        if(img(x_px(i),y_px(i)) == search_val)
            res = true;
            return
        end
    end

end

