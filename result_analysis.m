data_base = AuTPs;
data_subject = AuTPs;
class_1 = 2;
class_2 = 1;
class_also_include = 0;
color_1 = colors(2, :);
color_2 = colors(6, :);
% particles_f6d = particles_f6(ismember(step_results.step2.sub_classes{2}.classes, [1 2]));
% features_f6d = step_results.step2.sub_classes{2}.data(ismember(step_results.step2.sub_classes{2}.classes, [1 2]), :);
% classes_f6d = step_results.step2.sub_classes{2}.classes(ismember(step_results.step2.sub_classes{2}.classes, [1 2]));

features_base = data_base.features;
mean_norm = mean(features_base);
max_norm = max(features_base-mean(features_base));
features_base = (features_base-mean_norm)./max_norm;
classes_base = data_base.classification;
[mu, sigma, ~] = compute_distribution(features_base, classes_base);

features_subject = data_subject.features;
features_subject = (features_subject-mean_norm)./max_norm;
classes_subject = assignlabels(features_subject, mu, sigma);

%%
classes_f = classes_subject;
features_f = features_subject;
particles_f = data_subject.particle_shapes;
data_name = data_subject.data_name;

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

%%

classes_idx_subject = ismember(classes_subject, [class_1 class_2 class_also_include]);
particles_subject = data_subject.particle_shapes(classes_idx_subject);

%[mu, sigma] = compute_distribution(features_f6d_sub, classes_f6d_sub);
prob_1 = mvnpdf(features_subject(classes_idx_subject, :), mu(class_1, :), diag(sigma(class_1, :)));
prob_2 = mvnpdf(features_subject(classes_idx_subject, :), mu(class_2, :), diag(sigma(class_2, :)));
prob_1 = prob_1/max(prob_1);
prob_2 = prob_2/max(prob_2);
frac_1 = prob_1./(prob_1+prob_2);
frac_2 = prob_2./(prob_1+prob_2);

sample_ratio = 10;
colored_particles_f6d = cell(ceil(length(particles_subject)/sample_ratio), 1);
j = 1;

for i = 1:sample_ratio:length(particles_subject)
    
    colored_particles_f6d{j} = cat(3,...
        particles_subject{i}*(frac_1(i)*color_1(1)+frac_2(i)*color_2(1))...
        , particles_subject{i}*(frac_1(i)*color_1(2)+frac_2(i)*color_2(2))...
        , particles_subject{i}*(frac_1(i)*color_1(3)+frac_2(i)*color_2(3)));
    j = j+1;
    
end

[~, idx] = sort(frac_2(1:sample_ratio:length(particles_subject)));

figure;montage(colored_particles_f6d(idx), 'ThumbnailSize', [1 1]*64);
%%
particles_subject = data_subject.particle_shapes;
prob_1 = mvnpdf(features_subject, mu(class_1, :), diag(sigma(class_1, :)));
prob_1 = prob_1/max(prob_1);
prob_2 = 1-prob_1;
figure;
histogram(prob_1,50,'Normalization','probability')
sample_ratio = 10;
colored_particles_f6d = cell(ceil(length(particles_subject)/sample_ratio), 1);
j = 1;

for i = 1:sample_ratio:length(particles_subject)
    
    colored_particles_f6d{j} = cat(3,...
        particles_subject{i}*(prob_1(i)*color_1(1)+prob_2(i)*color_2(1))...
        , particles_subject{i}*(prob_1(i)*color_1(2)+prob_2(i)*color_2(2))...
        , particles_subject{i}*(prob_1(i)*color_1(3)+prob_2(i)*color_2(3)));
    j = j+1;
    
end

[~, idx] = sort(prob_2(1:sample_ratio:length(particles_subject)));

figure;montage(colored_particles_f6d(idx), 'ThumbnailSize', [1 1]*64);

%%
features_toplot = {features_AuTPs, features_H2O2, features_HNO3};
feature_names = {'Area', 'Eccentricity', 'Aspect Ratio', 'Circularity'};

for i=1:4
    maxs = [0 0 0];
    mins = [0 0 0];
    for n=1:3
        maxs(n) = max(features_toplot{j}(:,i));
        mins(n) = min(features_toplot{j}(:,i));
    end
    max_ite = max(maxs);
    min_ite = max(mins);
    figure
    hold on
    for j=1:3
        features_ite = features_toplot{j}(:,i);
        h = histogram(features_ite, 'BinWidth', 0.04*(max_ite-min_ite), 'FaceColor', colors(j, :), 'EdgeAlpha', 0.25,'Normalization','probability');
    end
    hold off
    title(feature_names{i})
    legend({'No etchant', 'H2O2', 'HNO3'})
end


%%
datasets = {AuTPs, AuTPs_10H2O2, AuTPs_50H2O2, AuTPs_100H2O2};
H2O2 = [0 10 50 100];
features_pool = [];
for i = 1:length(datasets)
    features_pool_ite = [datasets{i}.features H2O2(i)*ones([length(datasets{i}.features), 1])];
    features_pool = [features_pool; features_pool_ite];
end

features_pool_norm = (features_pool-mean(features_pool))./max(features_pool-mean(features_pool));
[results, step_results] = iterativeclustering(features_pool_norm, 5);
%%
particles_pool = [];
for i = 1:length(datasets)
    particles_pool = [particles_pool datasets{i}.particle_shapes];
end