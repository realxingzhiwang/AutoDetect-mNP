%%
%Initializing
addpath(genpath('Hu Moments'), 'iterativeclustering',genpath('loadEMimages'),'naivebayes','ruecs')
%%
%Find all images from a folder, make random split if desired
folder = uigetdir;
directs = dir(fullfile(folder, '*.dm4'));
%directs = dir([folder '\*.png']); %Defualt input format is .dm4
names = {directs.name};
%%
%Load images can calculated features
image = {};
image_bw = {};                                                                                                                                                                                                              
features_org = [];
particles = {};
particles_list = {};
% particles_org = {};
% centroids = [];
% orientation = {};
moments = [];
codes = [];
image_loading = @ReadDMFile; %@ReadDMFile for dm4, @loadtiff for other formats
image_segmentation = @imagekmeans; %@imagekmeans for mNPs, @combinedthresh for QDs, @identity if inputs are binary images
[~, ~, ~, ~, unit] = loadEMimages(fullfile(folder,names{1}), image_loading, image_segmentation);
for i = 1:length(names)
    [image{i}, image_bw{i}, features_ite, particles_ite, unit, particles_org_ite, moments_ite, centroids_ite, orientation_ite]...
        = loadEMimages(fullfile(folder,names{i}), image_loading, image_segmentation); 
    features_org = [features_org; features_ite];
    particles = [particles particles_ite{1}];
    particles_list = [particles_list particles_ite{2}];
%     particles_org = [particles_org particles_org_ite];
    moments = [moments; moments_ite];
%     centroids = [centroids; centroids_ite];
%     orientation = [orientation orientation_ite];
    codes = [codes; i*ones(length(particles_org_ite), 1)];
end

%Plot particles at there original positions (for making overlay images)
%particles_full = cell(size(particles_list));
%for i = 1:length(particles_list)
%    particles_full{i} = masking(image_bw{codes(i)}, particles_list{i});
%end

%%
%Filtering non-convex particles
nonoverlapping =...based on solidity and convexity
    features_org(:,end)>0.9 & features_org(:,7)>0.95;%0.9/0.95
overlapping = ~nonoverlapping;

features_nonoverlapping = features_org(nonoverlapping, :);
features_nonoverlapping = [features_nonoverlapping(:,[1 2 6]) 1./(moments(nonoverlapping,1)*2*pi)]; %area, eccentricity, aspect ratio, circularity
particles_plot = particles(nonoverlapping);
particles_list_plot = particles_list(nonoverlapping);

%%
%Perform classification
features_norm = (features_nonoverlapping-mean(features_nonoverlapping))./max(features_nonoverlapping-mean(features_nonoverlapping));
[results, step_results] = iterativeclustering(features_norm, 5);

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
overlapping_codes = codes(overlapping);
overlapping_list = particles_list(overlapping);
resolved_codes = [];
resolved_list = cell([]);

parfor i = 1:N
   markers{i} = ruecs(particles_ol{i}, (10/unit{1})^2);
   [markers_dil{i}, overlay{i}] = dilmarkers(markers{i}, particles_ol{i});
   if ~isempty(markers_dil{i})
       resolved_codes = [resolved_codes; overlapping_codes(i)*ones(length(markers_dil{i}), 1)];
       resolved_list = [resolved_list repmat(overlapping_list(i), 1, length(markers_dil{i}))];
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
codes_all = [codes(nonoverlapping); resolved_codes];
particles_for_overlay = [particles_plot resolved_markers];
particles_list_all = [particles_list_plot resolved_list];

%% Summarizing output data
summary = struct;
summary.particle_shapes = particles_all;
summary.features = features_all; %In the order: Area, Eccentricity, Aspect Ratio, Circularity
summary.classification = class_idx_all;
summary.original_images = image;
summary.binarized_images = image_bw;
summary.codes = codes_all;
summary.data_name = names{1};

%% Results visualization
colors = [68 133 255;
        133 255 68;
        255 68 133;
        255 157 37;
        157 37 255;
        37 255 157;
        176 176 176
        117 138 155]/255*0.9;
    
handles = cell(1, 4);
feature_names = {'Area', 'Eccentricity', 'Aspect Ratio', 'Circularity'};
summary_f = summary;
classes_f = summary_f.classification;
features_f = summary_f.features;
particles_f = summary_f.particle_shapes;
data_name = summary_f.data_name;

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
    title(data_name)
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
    title([data_name ' Class ' num2str(j)])
end

classes_counts = zeros(1, max(classes_f));
classes_labels = cell(1, max(classes_f));
for k = 1:max(classes_f)
    classes_counts(k) = sum(classes_f==k);
    classes_labels{k} = ['Class ' num2str(k)];
end

figure
ax = gca();
pie(ax, classes_counts)
legend(classes_labels)
title(data_name)
ax.Colormap = colors(1:max(classes_f), :);

% figure
% montage(colored_particles_all, 'ThumbnailSize', [1 1]*256)

%%
%image_w_class = cell(size(summary.binarized_images));
image_w_class_dir = 'D:\Wang_data\image_w_class\';
classes_for_overlay = summary.classification;
parfor i = 1:length(summary.original_images)

    %R_layer = zeros(size(particles_for_overlay{summary.codes==i}));
    %G_layer = zeros(size(particles_for_overlay{summary.codes==i}));
    %B_layer = zeros(size(particles_for_overlay{summary.codes==i}));
    R_layer = zeros(size(summary.original_images{i}));
    G_layer = zeros(size(summary.original_images{i}));
    B_layer = zeros(size(summary.original_images{i}));
    
    for j = 1:max(classes_for_overlay)
        particles_overlay_ite = particles_for_overlay(summary.codes==i & classes_for_overlay==j);
        particles_list_overlay_ite = particles_list_all(summary.codes==i & classes_for_overlay==j);
        if ~isempty(particles_overlay_ite)
            %R_layer = R_layer + colors(j, 1)*sum(cat(3, particles_overlay_ite), 3);
            %G_layer = G_layer + colors(j, 2)*sum(cat(3, particles_overlay_ite), 3);
            %B_layer = B_layer + colors(j, 3)*sum(cat(3, particles_overlay_ite), 3);
            for n = 1:length(particles_overlay_ite)
                box_size = size(particles_overlay_ite{n});
                y_range = round(particles_list_overlay_ite{n}(1)):round(particles_list_overlay_ite{n}(1))+box_size(2)-1;
                x_range = round(particles_list_overlay_ite{n}(2)):round(particles_list_overlay_ite{n}(2))+box_size(1)-1;

                R_layer(x_range, y_range) = R_layer(x_range, y_range) + colors(j, 1)*particles_overlay_ite{n};
                G_layer(x_range, y_range) = G_layer(x_range, y_range) + colors(j, 2)*particles_overlay_ite{n};
                B_layer(x_range, y_range) = B_layer(x_range, y_range) + colors(j, 3)*particles_overlay_ite{n};
            end

        end
    end
    
    image_w_class = summary.original_images{i}*0.75 + uint8((cat(3, R_layer, G_layer, B_layer))*255);
    imwrite(image_w_class, [image_w_class_dir names{i} '.png']);
end

%%
function [ masked_img ] = masking(img, region)
    mask = zeros(size(img));
    mask(region) = 1;
    masked_img = logical(img.*mask);
end

function out = getVarName(var)
    out = inputname(1);
end