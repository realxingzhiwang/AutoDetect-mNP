%%
%Initializing
addpath(genpath('Hu Moments'), 'iterativeclustering','loadEMimages','naivebayes','ruecs')
%%
%Find all images from a folder, make random split if desired
folder = uigetdir;
directs = dir([folder '\*.dm4']);
%directs = dir([folder '\*.tif']); Defualt input format is .dm4
names = {directs.name};
%%
%Load images can calculated features
image = {};
image_bw = {};                                                                                                                                                                                                              
features_org = [];
particles = {};
% particles_list = {};
% particles_org = {};
% centroids = [];
% orientation = {};
moments = [];
codes = [];
[~, ~, ~, ~, unit] = loadEMimages(fullfile(folder,names{1}));
for i = 1:length(names)
    [image{i}, image_bw{i}, features_ite, particles_ite, unit, particles_org_ite, moments_ite, centroids_ite, orientation_ite]...
        = loadEMimages(fullfile(folder,names{i})); %Add '@loadtiff' as the second input if loading non .dm4 format is desired
    features_org = [features_org; features_ite];
    particles = [particles particles_ite{1}];
%     particles_list = [particles_list particles_ite{2}];
%     particles_org = [particles_org particles_org_ite];
    moments = [moments; moments_ite];
%     centroids = [centroids; centroids_ite];
%     orientation = [orientation orientation_ite];
    codes = [codes; i*ones(length(particles_org_ite), 1)];
end

%%
%Filtering non-convex particles
nonoverlapping =...based on solidity and convexity
    features_org(:,end)>0.9 & features_org(:,7)>0.95;%0.8/0.9
overlapping = ~nonoverlapping;

features_nonoverlapping = features_org(nonoverlapping, :);
features_nonoverlapping = [features_nonoverlapping(:,[1 2 6]) 1./(moments(nonoverlapping,1)*2*pi)]; %area, eccentricity, aspect ratio, circularity
particles_plot = particles(nonoverlapping);

%%
%Perform classification
features_norm = (features_nonoverlapping-mean(features_nonoverlapping))./max(features_nonoverlapping-mean(features_nonoverlapping));
[results, step_results] = iterativeclustering(features_norm, 5);

%%
%rUECS
particles_ol = particles(overlapping);
N = length(particles_ol);
Img = cell(N, 1);
markers = cell(N, 1);
cnt = zeros(N, 1);
overlay = cell(size(markers));
layers = cell(size(markers));
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

resolved_features_norm = (resolved_features-mean(features_nonoverlapping))./...
    max(resolved_features-mean(features_nonoverlapping));

class_idx_max = step_results.step1.classes;
[mu_max, sigma_max] = compute_distribution(step_results.step1.data, class_idx_max);

class_idx_resolved = assignlabels(resolved_features_norm,...
    mu_max, sigma_max);

particles_all = [particles_plot resolved_particles'];
features_all = [features_nonoverlapping; resolved_features];
class_idx_all = [class_idx_max; class_idx_resolved];
codes_all = [codes; resolved_codes];

%% Summarizing output data
summary = struct;
summary.particle_shapes = particles_all;
summary.features = features_all; %In the order: Area, Eccentricity, Aspect Ratio, Circularity
summary.classification = class_idx_all;
summary.original_images = image;
summary.binarized_images = image_bw;

%% Results visualization
colors = [68 133 255;
        133 255 68;
        255 68 133;
        255 157 37;
        157 37 255;
        37 255 157]/255*0.9;
    
handles = cell(1, 4);
feature_names = {'Area', 'Eccentricity', 'Aspect Ratio', 'Circularity'};
classes_f = summary.classification;
features_f = summary.features;
particles_f = summary.particle_shapes;

for p = 1:4
    
    figure
    hold on
    for j=1:max(classes_f)
        features_f_ite = features_f(classes_f==j,p);
        h = histogram(features_f_ite, 'BinWidth', 0.04*(max(features_f(:,p))-min(features_f(:,p))), 'FaceColor', colors(j, :), 'EdgeAlpha', 0.25);
        ax = axis;
        gau_curve = normpdf(linspace(min(features_f(:,p)), max(features_f(:,p)), 100),...
            mean(features_f_ite),...
            std(features_f_ite));
        plot(linspace(min(features_f(:,p)), max(features_f(:,p)), 100),...
            gau_curve/max(gau_curve)*max(h.Values), 'Color', colors(j, :)*0.8, 'LineWidth', 0.75);

    end
    hold off
    axis tight
    ytickformat('%.2f')
    yt = get(gca, 'ytick')/size(classes_f, 1);
    set(gca, 'FontSize', 22, 'YTickLabel', compose('%.2f', string(yt)))
    xlabel(feature_names{p}, 'FontSize', 22)
    ylabel('Normalized Counts', 'FontSize', 22)
    box on
    pbaspect([1 1 1])
    handles{p} = gca;
end

colored_particles_all = {};
for j = 1:max(classes_f)
    particles_plot_f = particles_f(classes_f==j);
    colored_particles_f = cell(size(particles_plot_f));
    for i = 1:length(colored_particles_f)
        colored_particles_f{i} = cat(3, particles_plot_f{i}*colors(j, 1),...
            particles_plot_f{i}*colors(j, 2), particles_plot_f{i}*colors(j, 3))*1.5;
    end
    figure
    if length(colored_particles_f)>=16 %Plotting 16 particles in each class for visualization
        montage(colored_particles_f(1:16), 'BorderSize', [1 1], 'ThumbnailSize', [1 1]*256)
    else
        montage(colored_particles_f(1:end), 'BorderSize', [1 1], 'ThumbnailSize', [1 1]*256)
    end
    colored_particles_all = [colored_particles_all colored_particles_f];
end

% figure
% montage(colored_particles_all, 'ThumbnailSize', [1 1]*256)