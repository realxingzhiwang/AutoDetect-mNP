function [dtimage] = fanalyzeimages(sourceDir, blackbackground)
%UNTITLED16 Summary of this function goes here
%   Detailed explanation goes here

if not(nargin);
    sourceDir = uigetdir(path, 'select source folder'); 
end 
%% get the file list from the directory%%
% sourceFiles = dir(fullfile(sourceDir, '*.png'));
% type = 'png'; 
% if isempty(sourceFiles)
%     sourceFiles = dir(fullfile(sourceDir, '*.dm3'));
%     type = 'dm3'; 
% end
% if isempty(sourceFiles) 
%     sourceFiles = dir(fullfile(sourceDir, '*.dm4')); 
%     type = 'dm4'; 
% end 
% if isempty(sourceFiles)
%     sourceFiles = dir(fullfile(sourceDir, '*.tif'));
%     type = 'tif'; 
% end

%trying again with different types of files in one file...but it would be
%better if we just didn't do that kind of shit (Chelsea)

sourceFiles1 = dir(fullfile(sourceDir, '*.png')); 
sourceFiles2 = dir(fullfile(sourceDir, '*.dm3'));
sourceFiles3 = dir(fullfile(sourceDir, '*.dm4'));
sourceFiles4 = dir(fullfile(sourceDir, '*.tif'));

%create empty struct
sourceFiles.name = []; 
sourceFiles.folder = []; 
sourceFiles.date = []; 
sourceFiles.bytes = []; 
sourceFiles.isdir = []; 
sourceFiles.datenum = []; 
sourceFiles.type = []; 

currentsize = 1; 
f = fieldnames(sourceFiles);

if ~isempty(sourceFiles1)
    addit = length(sourceFiles1); 
    
    for h = 1:addit
        sourceFiles1(h).type = 'png';  
    end 
    
    f = fieldnames(sourceFiles); 
    numims = length(sourceFiles1); 
    for j = 1:length(sourceFiles1)
        for i = 1:length(f) 
            sourceFiles(currentsize+j).(f{i}) = sourceFiles1(j).(f{i}); 
        end 
    end 
    
    currentsize = currentsize + addit;
end 

if ~isempty(sourceFiles2) 
    
    addit = length(sourceFiles2); 
    
    for h = 1:addit
        sourceFiles2(h).type = 'dm3'; 
    end 

   
    for j = 1:length(sourceFiles2)
        for i = 1:length(f) 
            sourceFiles(currentsize+j).(f{i}) = sourceFiles2(j).(f{i}); 
        end 
    end 
    
    currentsize = currentsize + addit;
end 

if ~isempty(sourceFiles3)
    addit = length(sourceFiles3); 
    %sourceFiles3(1:addit).type = 'dm4';  
    
    for h = 1:addit
        sourceFiles3(h).type = 'dm4'; 
    end 
   
    for j = 1:length(sourceFiles3)
        for i = 1:length(f) 
            sourceFiles(currentsize+j).(f{i}) = sourceFiles3(j).(f{i}); 
        end 
    end 
    
    currentsize = currentsize + addit;
end 

if ~isempty(sourceFiles4) 
    addit = length(sourceFiles4); 
    sourceFiles4(1:addit).type = 'tif';  
   
    for j = 1:length(sourceFiles4)
        for i = 1:length(f) 
            sourceFiles(currentsize+j).(f{i}) = sourceFiles4(j).(f{i}); 
        end 
    end 
    
end 

%taking out the arbitrary first entry
sourceFiles(1) = []; 

%%For each image file analyze the image%%
for j=1:length(sourceFiles)
               tic; 
%              try
                Iname = sourceFiles(j).name; 
                type = sourceFiles(j).type; 
                
                if (strcmp(type,'tif'))|(strcmp(type,'png')); 
                I = imread(fullfile(sourceDir, sourceFiles(j).name));
                I = rgb2gray(I);
                mag = getmag(Iname); 
                end 
                
                if strcmp(type,'dm3'); 
                dm3struc = DM3Import(fullfile(sourceDir, sourceFiles(j).name)); 
                Idouble = mat2gray(dm3struc.image_data', [dm3struc.intensity.lowlim(1), dm3struc.intensity.highlim(1)]);
                I = im2uint8(Idouble); 
                mag = 1/(dm3struc.xaxis.scale);        
                end 
                
                if strcmp(type,'dm4'); 
                [dm4array,dm4mag,~] = ReadDMFile(fullfile(sourceDir, sourceFiles(j).name)); 
                Idouble = mat2gray(dm4array); 
                I = im2uint8(Idouble); 
                mag = 1/(dm4mag); 
                end 
                
                dtimage(j).name = Iname; 
                dtimage(j).original = I; 
                dtimage(j).magnification = mag; 
               
                if blackbackground
                    [BW, filtered, thresh] = fthreshold1(imcomplement(I),0,1); %fthreshold1(image,auto (when false), filter (when true))
                else
                    [BW, filtered, thresh] = fthreshold1(I,0,1);
                end 
                 
                dtimage(j).edited = filtered; 
                dtimage(j).binary = BW; 
                dtimage(j).threshold = thresh; 
                
                jthImageProperties = fCalculateProperties(BW, mag); 
                
                dtimage(j).properties = jthImageProperties; 
                elapsed = toc; 
                fprintf('done with %f / %f images  time: %f seconds \n', j,length(sourceFiles),elapsed); 
                
%             catch
%                 disp(strcat('Problem with: ', sourceFiles(j).name));
%             end
end 
end

function [mag] = getmag(name) 
    %input; name of the sample as a string  
    %must be of format name$mag
    try
    C = strsplit(name,'$'); 
    mag = str2num(C{2}(1:(end-4))); 
    catch
        disp('ERROR: if using a non DM3 file, \n please append the image name with the pixel/nm ratio after a  $ \n samplename$ratio.png \n'); 
    end 
end
