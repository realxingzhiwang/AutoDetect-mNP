function [H] = bulkcrossentropy(data,labels)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

classes = unique(labels);

mu_bulk = mean(data);
sigma_bulk = sqrt(var(data));
[mu_classes, sigma_classes] = compute_distribution(data, labels);

H = zeros(length(classes), 1);

for i = 1:length(classes)
    H(i) = -mvnpdf(data(labels==classes(i), :), mu_classes(i, :), diag(sigma_classes(i, :)))'...
        *log(mvnpdf(data(labels==classes(i), :), mu_bulk, diag(sigma_bulk)))...
        /sum(labels==classes(i));
end

end

