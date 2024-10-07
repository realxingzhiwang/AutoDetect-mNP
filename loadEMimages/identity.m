function [BW_fill_filter] = identity(image_8bit)
%For testing only
%   Detailed explanation goes here
image_8bit = imclearborder(image_8bit); %Discard any partciles truncated by the border

BW_fill_filter = image_8bit~=0;

BW_fill_filter = bwmorph(BW_fill_filter,'close');
BW_fill_filter = bwmorph(BW_fill_filter,'open');

end
