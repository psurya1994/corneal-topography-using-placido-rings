function contour_on_dot = is_contour_on_dot(r_contour, r_matrix)

    contour_on_dot = false;
    
    for i = 1:length(r_contour)
        if(isnan(r_contour(i)))
            continue
        elseif(r_matrix(i, r_contour(i)) == 1)
            contour_on_dot = true;
            return
        end
    end
end