function [image_16bit, scale, unit] = loadtiff(file)
%For testing only
%   Detailed explanation goes here
image_16bit = imread(file);
scale = 0.1580; %Scale for mixture sample (17kx)
unit = 'nm';
if size(image_16bit, 3)>1
    image_16bit = rgb2gray(image_16bit);
end
end

