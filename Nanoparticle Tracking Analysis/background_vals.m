function [avg_gray_with_black] = background_vals(image)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
avg_gray_with_black = mean(image, 'all');
%avg_gray_no_black = mean(image(image > 0), 'all');
end