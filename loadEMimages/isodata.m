function [BW_fill_filter] = isodata(image_8bit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
image_I = 255-image_8bit; %Invert color
ContI = imclearborder(image_I); %Discard any partciles truncated by the border
%[BW, thresh] = edge(ContI, 'log', 0);
BW = imsegisodata(ContI); %Use k-means to segment image

%BW = logical(L-1);
%BW = L-1;
%BW = logical(1-BW);

BW_fill_filter = imfill(BW,4,'holes');
if any(BW_fill_filter(:)) %Filter based on particle sizes
    BW_fill_filter = bwareafilt(BW_fill_filter, [500 5000000]);
end    
BW_fill_filter = bwmorph(BW_fill_filter,'spur');
BW_fill_filter = bwmorph(BW_fill_filter,'majority');
BW_fill_filter = bwmorph(BW_fill_filter,'close');
BW_fill_filter = bwmorph(BW_fill_filter,'bridge');
BW_fill_filter = bwmorph(BW_fill_filter,'open');
BW_fill_filter = imfill(BW_fill_filter,4,'holes');
BW_fill_filter = imfilter(BW_fill_filter,fspecial('gaussian', [10 10]));
end