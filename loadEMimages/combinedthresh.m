function [BW_fill_filter] = combinedthresh(image_8bit)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

BW_fill_filter = fthreshold1(image_8bit,0,1);

end

