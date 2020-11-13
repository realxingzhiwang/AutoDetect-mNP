function [image_8bit, varargout] = loadEMimages(varargin)
%[image_8bit, image_binary, features, particles, scale, particles_org, moments] = loadEMimages(file (, loading_function))
%   loadEMimages loads TEM images from the directory specified by file
%   using function specified by loading_function. file is a string of the
%   directory of the TEM image file. loading_function is a function handle
%   specifying a function used to load the image files. By default,
%   loading_function is @ReadDMFile, which loads .dm3 and .dm4 formats. A
%   valid loading function must provide output in the form of [image, pixel
%   scale, unit], where iamge is an int matrix, pixel scale is a double,
%   and unit is a string (eg. 'nm'). loadEMimages outputs image_8bit,
%   image_binary, features, particles, scale, particles_org, moments, where
%   image_8bit is an int8 matrix representing original image in greyscale,
%   image_binary is a logical matrix representing the image segmented by
%   k-means clustering, with 1's being foreground and 0's being background.
%   features is a N-by-9 matrix, where N is the number of isolatable
%   particles in the segmented image. Each row of features represents a set
%   of 9 features computed by regionprops (specified below) for a specific
%   isolated particle. particles is a N-by-1 cell, with each element being
%   a logical matrix representing an isolated particle. unit is a 1-by-2
%   cell, with unit{1} being the pixel scale of the image, as a double, and
%   unit{2} being the unit os the scale, as a string. particles_org is a
%   N-by-1 cell, with each element being an int8 matrix representing the
%   cropped greyscale image containing one isolated particle. moments is a
%   N-by-7 matrix with each row being the 7 Hu's moment invariant
%   calculated for an isolated partcile.

if nargin==2
    file = varargin{1};
    loadimage = varargin{2};
elseif nargin==1
    file = varargin{1};
    loadimage = @ReadDMFile;
end


if nargin<3
    [image_32bit, scale, unit] = loadimage(file);
else
    image_32bit = varargin{1};
    scale = varargin{2};
    unit = varargin{3};
end


area_threshold = scale^2*0.5e4/10;
%area_threshold = 100;
%image_16bit = imadjust(image_16bit);
% image_8bit = imread(file);
% image_8bit = uint8(conv2(image_8bit(1:1536,1:1536), 1/25*ones(5)));
% scale = 1;
% unit = 'nm';

% if unit ~= 'nm'
%     disp('Unit not in nm');
% end

%image_16bit = imfilter(image_16bit,fspecial('gaussian'));
if length(size(image_32bit))>2 %Conversion from rgb to greyscale
    image_32bit = rgb2gray(image_32bit);
end
% if max(image_32bit(:))>257 %Conversion from int16 to int8
%     image_8bit = uint8(conv2(image_32bit/512, 1/25*ones(5))); %The factor 512 can be tuned based on the dataset
% else
%     image_8bit = uint8(image_32bit);
% end
image_8bit = im2uint8(mat2gray(image_32bit));
image_8bit = imadjust(image_8bit); %Adjust contrast

% image_dim = floor(size(image_8bit)*0.65);
% image_sections = uint8(zeros([image_dim, 4]));
% BW_sections = zeros([image_dim, 4]);
% BW_fill_filter = zeros(size(image_8bit));
% section_boundaries = {1:image_dim(1), 1:image_dim(2);...
%     size(image_8bit, 1)-image_dim(1)+1:size(image_8bit, 1),...
%     size(image_8bit, 2)-image_dim(2)+1:size(image_8bit, 2)};
% for i = 1:4
%    image_sections(:, :, i) = uint8(image_8bit(section_boundaries{1+(i>2), 1}, section_boundaries{2-mod(i, 2), 2}));
%    BW_sections(:, :, i) = imagekmeans(image_sections(:, :, i));
%    BW_fill_filter(section_boundaries{1+(i>2), 1}, section_boundaries{2-mod(i, 2), 2}) = ...
%        or(BW_fill_filter(section_boundaries{1+(i>2), 1}, section_boundaries{2-mod(i, 2), 2}),...
%        BW_sections(:, :, i));
% end

BW_fill_filter = imagekmeans(image_8bit);
%BW_fill_filter = combinedthresh(image_8bit);
%BW_fill_filter = bwmorph(BW_fill_filter,'hbreak');
%BW_fill_filter = ~BW_fill_filter;


% features
%area = regionprops(BW_fill_filter, 'Area');
%bb = regionprops(BW_fill_filter, 'BoundingBox');
%eccentricity = regionprops(BW_fill_filter, 'Eccentricity');
%extent = regionprops(BW_fill_filter, 'Extent');
%majoraxis = regionprops(BW_fill_filter, 'MajorAxisLength');
%minoraxis = regionprops(BW_fill_filter, 'MinorAxisLength');
%solidity = regionprops(BW_fill_filter, 'Solidity');
%perimeter = regionprops(BW_fill_filter, 'Perimeter');

% IDs
%centroid = regionprops(BW_fill_filter, 'Centroid');
%isolatedimage = regionprops(BW_fill_filter, 'Image');
%particles = {isolatedimage.Image};

regions = regionprops(BW_fill_filter, 'Area', 'BoundingBox',...
    'Eccentricity', 'Extent', 'MajorAxisLength', 'MinorAxisLength',...
    'Solidity', 'Perimeter', 'Centroid', 'Image', 'ConvexHull',...
    'Orientation', 'PixelIdxList'); %Compute features

particles = {regions.Image};
particles_full = {regions.PixelIdxList};
centroids_temp = [regions.Centroid];
centroids = reshape(centroids_temp, 2, length(centroids_temp)/2)';


if isempty(particles)
    
    if nargout ~= 0
    varargout{1} = zeros(size(image_8bit));
    end

    if nargout > 1
            varargout{2} = [];
    end

    if nargout > 2
            varargout{3} = {particles, {}};
    end

    if nargout > 3
            varargout{4} = {scale, unit};
    end
    
    if nargout > 4
            varargout(5:end) = cell(size(varargout(5:end)));
    end
    
    return
    
end

box = {regions.BoundingBox};
orientation = {regions.Orientation};
convex = {regions.ConvexHull};
convexperi = zeros(length(particles), 1);
particles_org = cell([1, length(particles)]);
extent = zeros(length(particles), 1);
moments = zeros(length(particles), 7);

warning('off','MATLAB:polyshape:repairedBySimplify')

for i = 1:length(particles)
    
    convexperi(i) = perimeter(polyshape(convex{i}));
    
    particles_org{i} = image_8bit...
        (ceil(box{i}(2)):ceil(box{i}(2))+box{i}(4)-1,...
        ceil(box{i}(1)):ceil(box{i}(1))+box{i}(3)-1);
    
    stats_rot = regionprops(imrotate(particles{i}, 90-orientation{i}),...
        'Extent');
    extent(i) = stats_rot.Extent;
    
    eta = SI_Moment(particles{i});
    moments(i, :) = Hu_Moments(eta);
    
end
    

%features = double([[area.Area]'*scale^2 ...
%    [eccentricity.Eccentricity]'...
%    [extent.Extent]'...
%    [majoraxis.MajorAxisLength]'*scale...
%    [minoraxis.MinorAxisLength]'*scale...
%    [majoraxis.MajorAxisLength]'./[minoraxis.MinorAxisLength]'...
%    [solidity.Solidity]'...
%    [area.Area]'./[perimeter.Perimeter]'.^2]);

features = double([[regions.Area]'*scale^2 ... area (nm^2)
    [regions.Eccentricity]'... eccentricity
    extent... extent
    [regions.MajorAxisLength]'*scale... major axis (nm)
    [regions.MinorAxisLength]'*scale... minor axis (nm)
    [regions.MajorAxisLength]'./[regions.MinorAxisLength]'... aspect ratio
    [regions.Solidity]'... solidity
    [regions.Area]'./[regions.Perimeter]'.^2 ... circularity
    convexperi./[regions.Perimeter]' ... convexity
    ]);


particles = particles(features(:,1)>area_threshold);
particles_full = particles_full(features(:,1)>area_threshold);
particles_org = particles_org(features(:,1)>area_threshold);
moments = moments(features(:,1)>area_threshold,:);
features = features(features(:,1)>area_threshold,:);
centroids = centroids(features(:,1)>area_threshold,:);


if nargout ~= 0
    varargout{1} = BW_fill_filter;
end

if nargout > 1
        varargout{2} = features;
end

if nargout > 2
        varargout{3} = {particles, particles_full};
end

if nargout > 3
        varargout{4} = {scale, unit};
end

if nargout > 4
    varargout{5} = particles_org;
end

if nargout > 5
    varargout{6} = moments;
end

if nargout > 6
    varargout{7} = centroids;
end

if nargout > 7
    varargout{8} = orientation;
end

end

