function [image_16bit, scale, unit] = loadtiff(file)
%For testing only
%   Detailed explanation goes here
image_16bit = imread(file);
scale = 0.1580; %Scale for mixture sample (17kx)
unit = 'nm';
end

