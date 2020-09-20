function [image_16bit, scale, unit] = loadtiff(file)
%For testing only
%   Detailed explanation goes here
image_16bit = imread(file);
scale = 1;
unit = 'nm';
end

