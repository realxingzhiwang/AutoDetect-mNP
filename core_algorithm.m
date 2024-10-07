function [summary] = core_algorithm(features_org,particles,moments,area_threshold, unit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Filtering non-convex particles
nonoverlapping =...based on solidity and convexity
    features_org(:,end)>0.9 & features_org(:,7)>0.95;%0.9/0.95
overlapping = ~nonoverlapping;

features_nonoverlapping = features_org(nonoverlapping, :);
features_nonoverlapping = [features_nonoverlapping(:,[1 2 6]) 1./(moments(nonoverlapping,1)*2*pi)]; %area, eccentricity, aspect ratio, circularity
particles_plot = particles(nonoverlapping);

%%
%Perform classification
features_norm = (features_nonoverlapping-mean(features_nonoverlapping))./max(features_nonoverlapping-mean(features_nonoverlapping));
[results, step_results] = iterativeclustering(features_norm, 2);

%%
%rUECS
%particles_ol = particles_full(overlapping);
particles_ol = particles(overlapping);
N = length(particles_ol);
Img = cell(N, 1);
markers = cell(N, 1);
cnt = zeros(N, 1);
overlay = cell(size(markers));
layers = cell(size(markers));
markers_dil = cell(size(markers));

parfor i = 1:N
   markers{i} = ruecs(particles_ol{i}, area_threshold);
   [markers_dil{i}, overlay{i}] = dilmarkers(markers{i}, particles_ol{i});
end

% 
% figure
% montage(overlay)

%%
%Computing features for resolved particles

scale_markers = unit{1};
resolved_markers = {};

parfor ite = 1:length(markers_dil)
    resolved_markers = [resolved_markers markers_dil{ite}];
end

resolved_features = zeros(length(resolved_markers), 4);
resolved_particles = cell(length(resolved_markers), 1);
resolved_axis = zeros(length(resolved_markers), 2);

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
    resolved_axis(ite, :) = [[props.MajorAxisLength],[props.MinorAxisLength]]*scale_markers;
end


%%
%Classifying resolved particles

resolved_features_norm = (resolved_features-mean(features_nonoverlapping))./...
    max(resolved_features-mean(features_nonoverlapping));

%class_idx_max = step_results.step1.classes;
class_idx_max = results.classes;
[mu_max, sigma_max] = compute_distribution(step_results.step1.data, class_idx_max);

class_idx_resolved = assignlabels(resolved_features_norm,...
    mu_max, sigma_max);

particles_all = [particles_plot resolved_particles'];
features_all = [features_nonoverlapping; resolved_features];
class_idx_all = [class_idx_max; class_idx_resolved];


%% Summarizing output data
summary = struct;
summary.particle_shapes = particles_all;
summary.features = features_all; %In the order: Area, Eccentricity, Aspect Ratio, Circularity
summary.classification = class_idx_all;
end