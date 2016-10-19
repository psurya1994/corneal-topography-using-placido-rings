function [ r ] = points2r(points, center)
%POINTS2R Summary of this function goes here
%   Detailed explanation goes here
    
    % r (segment, ring_number)
    r = zeros(length(points), length(points{1}));
    for i = 1:length(points) % segment
        vector = points{i}-repmat(center,length(points{i}),1);
        r(i,:) = sqrt(vector(:,1).^2 + vector(:,2).^2);
    end

end

