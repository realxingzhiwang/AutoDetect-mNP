function [markers] = ruecs(Img, area_threshold, cnt)
% Ultimate erosion of convex shapes adapted from Park, C. et al. Ding, Y.
%    IEEE Trans. Pattern Anal. Mach. Intell. 2013, 35, 669 - 680.
%   Recursive version

%ruecs Recursive ultimate erosion of convex shapes
%   [ markers ] = ruecs(img) returns all convex markers and their
%   corresponding number of itrations of erosions. To begin, input Img as a
%   logical matrix of an image of overlapping particles, and area_threshold
%   as the minimum plausible size of isolated particles as a double.
%   markers is a N-by-1 cell array, where N is the number of nonoverlapping
%   particles isolated from the input image. markers{n}.image is a
%   logical matrix representing the n-th isolated particle.
%   markers{n}.area is the area of the n-th isolated particle.
%   markers{n}.cnt is the number of iterations required to isolate that
%   particle.

if ~iscell(Img)
    s_init = regionprops(Img, 'Area');
    Img_tmp = Img;
    Img = struct;
    Img.image = Img_tmp;
    Img.init_area = s_init.Area;
    Img.area = s_init.Area;
    if nargin >= 3
        Img.cnt = cnt;
    else
        Img.cnt = 0;
    end
    Img.isconvex = false;
    Img.keep = true;
    Img = {Img};
end


se1    = strel('disk', 1 );  
se2    = strel('arbitrary', ones(2,2));
%se1 = se2;
sesize = 4; 


for main_ite = 1:length(Img)
    image = Img{main_ite}.image;
    A_0 = Img{main_ite}.init_area;
    warning('off','MATLAB:polyshape:repairedBySimplify')
    if ~Img{main_ite}.isconvex
        s = regionprops(image, 'PixelIdxList', 'ConvexArea', 'Area', 'ConvexHull', 'Perimeter');
        %cc = concavity(image);
        %if (max(cc) > 0.1 || 1 - s.Area / s.ConvexArea > 0.1)
        if (1 - s.Area / s.ConvexArea > 0.1) || (perimeter(polyshape(s.ConvexHull))/s.Perimeter < 0.9)%0.1/0.9 for TEM
                %&& s.MinorAxisLength > 9
            if mod(Img{main_ite}.cnt, 2) == 1
                se = se1;
            else
                se = se2;
            end
            Img_eroded = imerode(image, se);
            Img_eroded = imopen(Img_eroded, se1);
            s_eroded = regionprops(Img_eroded, 'PixelIdxList', 'Area');
            if isempty(s_eroded)
                Img{main_ite}.keep = false;
                continue
            end
            Imgs = {s_eroded.PixelIdxList};
            Areas = {s_eroded.Area};

            Img{main_ite}.image = masking(image, Imgs{1});
            Img{main_ite}.cnt = Img{main_ite}.cnt + 1;
            Img{main_ite}.area = Areas{1};
            if Areas{1} < 0.1*A_0 || Areas{1} < 25 %25 %0.3
                Img{main_ite}.keep = false;
            end
            
            if length({s_eroded.Area}) > 1  
                Img{main_ite}.area = Areas{1};
                A_0 = Areas{1};
                
                for sub_ite = 2:length(Imgs)
                    Img = [Img, ruecs(masking(image, Imgs{sub_ite}),...
                        area_threshold, Img{main_ite}.cnt)];
                end
            end

        else
            Img{main_ite}.isconvex = true;
        end
    end
        
end

isconvex = false(length(Img), 1);
keep = false(length(Img), 1);

for end_ite = 1:length(Img)
    if Img{end_ite}.area < area_threshold
        Img{end_ite}.keep = false;
    end
    isconvex(end_ite) = Img{end_ite}.isconvex;
    keep(end_ite) = Img{end_ite}.keep;
end

markers = Img(keep);

if isempty(markers)
    return
elseif sum(isconvex)<length(isconvex)
    markers = ruecs(markers, area_threshold);
end


function [ masked_img ] = masking(img, region)
    mask = zeros(size(img));
    mask(region) = 1;
    masked_img = logical(img.*mask);
end

end
