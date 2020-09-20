function [labels] = assignlabels(data, mu, sigma)
%[labels] = assignlabels(data, mu, sigma)
%   assignlabels assign a class labels for each datapoints in data, based
%   on the likelihood computed by assuming Gaussian class distributions
%   defined by mu and sigma. Inputs and outputs take the same format as
%   those specified in naivebayes.

likelihoods = zeros(size(data, 1), size(mu, 1));

for i = 1:size(mu, 1)
    if det(diag(sigma(i, :)))<=0
        likelihoods(:, i) = -Inf; %Eliminates empty classes and classes containing only one data point
    else
        likelihoods(:, i) = log(mvnpdf(data, mu(i, :), diag(sigma(i, :))));
    end
end

[~, labels] = max(likelihoods, [], 2);





end

