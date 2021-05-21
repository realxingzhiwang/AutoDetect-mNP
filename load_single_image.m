addpath(genpath('Hu Moments'), 'iterativeclustering',genpath('loadEMimages'),'naivebayes','ruecs')
[file_test, path_test] = uigetfile(('*.*'));
%[m, sx, units]=ReadDMFile(fullfile(path, file));
%imshow(uint8(m/256))
%%
[image_test, image_bw_test, features_test, particles_all_test, unit_test, particles_org_test, moments_test, centroids_test] = loadEMimages(fullfile(path_test, file_test),@loadtiff);
%image = imread(fullfile(path, file));
%image_pooled = averagepooling(image, [5 5]);
%image_conv = uint8(conv2(image, 1/25*ones(5)));
particles_test = [particles_all_test{1}];
figure
subplot(1,2,1)
imshow(image_test)
title('Original image')
subplot(1,2,2)
imshow(image_bw_test)
title('Binary image')