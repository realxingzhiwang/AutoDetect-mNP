function [mu, sigma, pi] = compute_distribution(data, labels)
%[mu, sigma, pi] = compute_distribution(data, labels)
%   compute_distribution computes mean (mu), standard deviation (sigma),
%   likelihood (pi) for each class, specified by labels, in data. The
%   format of inputs and outputs are the same as those specified in
%   naivebayes.

K = length(unique(labels));
mu = zeros(K, size(data, 2));
sigma = zeros(K, size(data, 2));
pi = zeros(K, 1);


for i = 1:K
    if ~isempty(data(labels==i, :))
        mu(i, :) = mean(data(labels==i, :));
        sigma(i, :) = sqrt(var(data(labels==i, :)));
        pi(i) = sum(labels==i)/length(labels);
    end
end



end

