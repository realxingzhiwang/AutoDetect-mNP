data_base = AuTPs_10min;
data_subject = AuTPs_10min;
class_1 = 2;
class_2 = 6;
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