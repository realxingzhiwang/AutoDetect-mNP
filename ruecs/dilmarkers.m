function [ markers_dil, overlay ] = dilmarkers(markers, image)
%dilmarkers Dilates markers obtained from ruecs
%   [ markers_dil, overlay ] = dilmarkers(markers) outputs a cell array containing
%   all the dilated markers. markers is a cell array output by ruecs.
%   markers_dil is a cell array with the same size as markers.
%   markers_dil{n} is a logical matrix representing the n-th resolved
%   particle. overlay is a 3D matrix of an rgb image showing the originial
%   marker, resolved markers and dilated markers overlaid together.
%   [ markers_dil, image_rec ] = dilmarkers(markers, image) reconstructs
%   images by overlaying dilated markers on the original image. This
%   functionality is planned but not yet implemented.

if isempty(markers)
    if nargout==1
        markers_dil = {};
    elseif nargout==2
        markers_dil = {};
        overlay = cat(3, zeros(size(image)), zeros(size(image)),...
            uint8(image)*255);
    end
    return
end


se1    = strel('disk', 1 );  
se2    = strel('arbitrary', ones(2,2));

markers_org = cell(size(markers));
markers_dil = cell(size(markers));

for i = 1:length(markers)
    markers_org{i} = markers{i}.image;
    markers_dil{i} = markers{i}.image;
    cnt = markers{i}.cnt;
    
    for j = cnt:-1:1
        if mod(j, 2) == 1
            se = se1;
        else
            se = se2;
        end
        markers_dil{i} = imdilate(markers_dil{i}, se);
    end

end

if nargin > 1
    
    outlines = false([size(image), length(markers)]);
    
    for i = 1:length(markers)
        outlines(:, :, i) = imdilate(bwmorph(markers_dil{i},'remove'),...
            strel('disk', 1));
    end
    
    outlines = sum(outlines, 3)>0;
    layer1 = uint8(sum(cat(3, markers_org{:}), 3)>0)*255;
    layer2 = uint8(sum(cat(3, markers_dil{:}), 3)>0)*255;
    layer3 = uint8(image)*255;
    
    layer1(outlines) = 255;
    layer2(outlines) = 0;
    layer3(outlines) = 128;
    
    overlay = cat(3, layer1, layer2, layer3);
end
    

end

