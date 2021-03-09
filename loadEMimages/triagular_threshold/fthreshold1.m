function [BWfinal,imagefiltered,level] = fthreshold1(image,manualthreshold,usefilter)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
    if (usefilter)
        imagefiltered = filter(image,5,7,true); 
    else
        imagefiltered = image; 
    end 

  
    if manualthreshold>0
        level = manualthreshold; 
    else 
        [level, EM] = graythresh(imagefiltered);     
        histogram = imhist(255-imagefiltered); 
        level = (1-triangle_th(histogram,256) + level)/2;  
    end 
    
% fprintf('threshold level is %f with effectiveness %f', level, EM); 

    BW = im2bw(imagefiltered, level); %optimum at 0.575
    BW = ~BW; 
    
    %%THIS LINE DETERMINES THE border CUT, to include the edge (e.g. for packing fraction calc) comment this line and uncomment the line below it
    BWcut = imclearborder(fwatershed(BW)); 
    % BWcut = fwatershed(BW); 

    BWfinal = repeatopen(BWcut,1,10,10); 
    % figure(5); 
    % imshow(BW);
%     plotboundaries(imagefiltered,BWfinal); 
%     plotboundaries(image,BWfinal);

% figure; 
% fprintf('threshold'); 
% imhist(imagefiltered); 
% figure;
% imshow(imagefiltered)
% figure; 
% imshow(BW);
% figure;
% Igs = grayslice(imagefiltered,16);
% imshow(Igs, jet(16)); 

end

function [filtered] = filter(I, n,s,use)
%Inputes:  I = image to be filtered, n = repetitions of filter s = kernal
%size in pixels
for i = 1:n
    if i == 1; 
          %I = medfilt2(I,[10,10]*3,'symmetric'); 
          I = imgaussfilt(I,12);
    end
    if i > 1; 
        %I = medfilt2(I,[s,s],'symmetric');
        I = imgaussfilt(I,s/4);
    end
end 
if (use)
se = strel('disk',100);
%I = imcomplement(imtophat(imcomplement(I),se));
filtered = adapthisteq(I); 
else
filtered = I; 
end
end

function [] = plotboundaries(I,BW)
%plots selected boundaries over original image
B = bwboundaries(BW); 
figure; 
 imshow(I); 
 hold on 
 for k=1:length(B)
            b = B{k};
            plot(b(:,2),b(:,1),'r','LineWidth',1);
 end
end 

function [cut] = fwatershed(BW)
    D = -bwdist(~BW);
    D = imhmin(D,10);
    % figure; 
    % imshow(D,[]); 
    L = watershed(D);
    BW(L == 0) = 0;
    cut = BW; 
end

function [Id] = repeatopen(I,j,e,d)
  sed = strel('disk', d);
  see = strel('disk',e);
  for i = 1:j; 
     I = imerode(I,see); 
     I = imdilate(I,sed);
  end 
  Id = I; 
end 