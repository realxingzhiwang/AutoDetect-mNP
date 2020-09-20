function [BW_fill_filter] = imagekmeans(image_8bit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
image_I = 255-image_8bit; %Invert color
ContI = imclearborder(image_I); %Discard any partciles truncated by the border
%[BW, thresh] = edge(ContI, 'log', 0);
L = imsegkmeans(ContI,2); %Use k-means to segment image
BW = false(size(L));
if mean(ContI(L==1))>mean(ContI(L==2)) %Set foreground to 1 and background to 0
    BW(L==1) = true;
else
    BW(L==2) = true;
end
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

