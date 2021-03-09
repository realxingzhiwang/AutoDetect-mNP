function [propertystruc] = fCalculateProperties(binaryimage,mag)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here
%area in pixels, perimeter
%mag = pixels/nm
C = regionprops(binaryimage, 'Area', 'Perimeter', 'MajorAxisLength','MinorAxisLength','Centroid','Orientation','Extrema', 'BoundingBox'); 

for j = 1:length(C)
C(j).nmArea = C(j).Area /(mag^2); 
C(j).diameter = 2*sqrt(C(j).nmArea/pi); 
%Calculate circularity 
%Circ = 4*pi*area/perimeter^2
C(j).Circ = C(j).Area /(C(j).Perimeter+3.14)^2*4*pi; 
%calculate roundess
%Round = 4*area/(pi*majoraxis^2)
C(j).Round = C(j).Area /(C(j).MajorAxisLength.^2)*4/pi; 

    %%%% For ellipse keep the code below comements and line 33 uncommented. For sphere do the reverse
    % distanceprev = 0; 
    % for i = 1:8
    %    for g = (i+1):8 
    %        x1 = C(j).Extrema(i,1); 
    %        x2 = C(j).Extrema(g,1); 
    %        y1 = C(j).Extrema(i,2); 
    %        y2 = C(j).Extrema(g,2); 
    %     distance = sqrt((x1-x2)^2 + (y1-y2)^2); 
    %         if(distance>distanceprev)
    %             distanceprev = distance; 
    %         end
    %     end    
    % end
    % C(j).nmFeret = distanceprev/mag; 
    C(j).nmFeret = C(j).MajorAxisLength/mag; 
end

propertystruc = C; 

end

