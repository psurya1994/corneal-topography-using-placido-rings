function [ mask ] = ellipse_mask(imageSize, center, a, b, phi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    [xx,yy] = ndgrid((1:imageSize(1))-center(1),(1:imageSize(2))-center(2));
    mask = ((xx.^2/a^2 + yy.^2/b^2)<1);
    
end

