function cen_mmt = Centr_Moment(image,mask,p,q)
% Function to calculate the central moment of interested image region
% Author:   Vishnu Muralidharan
% University of Alabama in Huntsville

% Inputs:   image: image: input image for which moments need to be calculated
%           mask: specifying this allows you to calculate moments for a
%           specified region
%           p,q: order of moments to be calculated
% Outputs:  cen_mmt = central moment of the specifed order fot the image
% Reference:  Visual Pattern Recognition by Moment Invariants


if ~exist('mask','var')
    mask = ones(size(image,1),size(image,2)); %if mask is not spcified, select the whole image
end

image = double(image);

%moments necessary to compute components of centroid
m10 = moment(image,mask,1,0); 
m01 = moment(image,mask,0,1);
m00 = moment(image,mask,0,0);

%components of centroid
x_cen = floor(m10/m00);
y_cen = floor(m01/m00);

cen_mmt =0;

a = size(mask, 1);
b = size(mask, 2);
i_mat = repmat([1:a]', 1, b) - x_cen;
j_mat = repmat([1:b], a, 1) - y_cen;
cen_mmt = sum(double(image .* (i_mat .^ p) .* (j_mat .^ q) .* mask), 'all');