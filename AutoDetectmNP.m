%%
%Initializing
addpath(genpath('Hu Moments'), 'iterativeclustering','loadEMimages','naivebayes','ruecs')
%%
%Find all images from a folder, make random split if desired
folder = uigetdir;
directs = dir([folder '\*.dm4']);
%directs = dir([folder '\*.tif']);
names = {directs.name};
%%
%Load images can calculated features
image = {};
image_bw = {};                                                                                                                                                                                                              
features_org = [];
particles = {};
particles_list = {};
% particles_full = {};
particles_org = {};
centroids = [];
orientation = {};
moments = [];
codes = [];
[~, ~, ~, ~, unit] = loadEMimages(fullfile(folder,names{1}));
for i = 1:length(names)%selectedims(1:20)
    [image{i}, image_bw{i}, features_ite, particles_ite, unit, particles_org_ite, moments_ite, centroids_ite, orientation_ite]...
        = loadEMimages(fullfile(folder,names{i}));
    features_org = [features_org; features_ite];
    particles = [particles particles_ite{1}];
    particles_list = [particles_list particles_ite{2}];
%     for j = 1:length(particles_ite{2})
%         particles_full = [particles_full masking(image_bw{i}, particles_ite{2}{j})];
%     end
    particles_org = [particles_org particles_org_ite];
    moments = [moments; moments_ite];
    centroids = [centroids; centroids_ite];
    orientation = [orientation orientation_ite];
    codes = [codes; i*ones(length(particles_org_ite), 1)];
end

%%
%Filtering non-convex particles
nonoverlapping =...based on solidity and convexity
    features_org(:,end)>0.9 & features_org(:,7)>0.95;
overlapping = ~nonoverlapping;

features = features_org(nonoverlapping, :);
features = [features(:,[1 2 6]) 1./(moments(nonoverlapping,1)*2*pi)]; %area, eccentricity, aspect ratio, circularity
particles_plot = particles(nonoverlapping);
features_plot = moments(nonoverlapping, :);

%%
%Perform classification
features_norm = (features-mean(features))./max(features-mean(features));
[results, step_results] = iterativeclustering(features_norm, 5);

%%
%rUECS
particles_ol = particles(overlapping);
%particles_ol = particles_plot(class_idx_em==2);
%particles_ol = all_ol;
N = length(particles_ol);
Img = cell(N, 1);
markers = cell(N, 1);
cnt = zeros(N, 1);
overlay = cell(size(markers));
layers = cell(size(markers));
%outlines_x = cell(size(markers));
%outlines_y = cell(size(markers));
markers_dil = cell(size(markers));
overlapping_codes = codes(overlapping);
resolved_codes = [];

parfor i = 1:N
   markers{i} = ruecs(particles_ol{i}, (10/unit{1})^2);
   [markers_dil{i}, overlay{i}] = dilmarkers(markers{i}, particles_ol{i});
   if ~isempty(markers_dil{i})
       resolved_codes = [resolved_codes; overlapping_codes(i)*ones(length(markers_dil{i}), 1)];
   end
end

figure
montage(overlay)

%%
%Computing features for resolved particles

scale_markers = unit{1};
resolved_markers = {};

parfor ite = 1:length(markers_dil)
    resolved_markers = [resolved_markers markers_dil{ite}];
end

resolved_features = zeros(length(resolved_markers), 4);
resolved_particles = cell(length(resolved_markers), 1);

parfor ite = 1:length(resolved_markers)
    props = regionprops(resolved_markers{ite},'Image', ...
        'Area', 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength');
    eta = SI_Moment(resolved_markers{ite});
    resolved_moments = Hu_Moments(eta);
    resolved_features(ite, :) = double([[props.Area]'*scale_markers^2, ...
        [props.Eccentricity]', ...
        [props.MajorAxisLength]'./[props.MinorAxisLength]', ...
        1./(resolved_moments(1)*2*pi)]);
    resolved_particles{ite} = props.Image;
end

%%
%Classifying resolved particles

resolved_features_norm = (resolved_features-mean(features))./...
    max(resolved_features-mean(features));

class_idx_max = class_idx_em(:, idx_max);
params_max = params_all{idx_max};

class_idx_resolved = assignlabels(resolved_features_norm,...
    params_max.mu, params_max.sigma);

for j=1:max(class_idx_max)
    figure('NumberTitle', 'off', 'Name', ['Class' num2str(j)]);
    montage(resolved_particles(class_idx_resolved==j), 'BorderSize', [1 1])
    figure('NumberTitle', 'off', 'Name', ['Class' num2str(j)]);
    for i = 1:size(features, 2)
        subplot(2,2,i);
        hist(resolved_features(class_idx_resolved==j,i),50);
    end
end

particles_all = [particles_plot resolved_particles'];
class_idx_all = [class_idx_max; class_idx_resolved];
codes_all = [codes; resolved_codes];