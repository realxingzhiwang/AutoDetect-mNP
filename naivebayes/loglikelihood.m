function [L] = loglikelihood(data, labels, mu, sigma, pi)
%[L] = loglikelihood(data, labels, mu, sigma, pi)
%   loglikelihood computes the total log likelihood of all classes,
%   specified by labels, in data, assuming Gaussian distribution defined by
%   mu and sigma. Inputs have the same format as those specified in
%   naivebayes. L is the total log likelihood in the format of a double.

L = 0;

for i = 1:length(pi)
    
    %if ~isempty(data(labels==i,:))
    if det(diag(sigma(i, :)))>0
        S = diag(sigma(i, :));
    
        P = mvnpdf(data(labels==i,:), mu(i, :), S);
        %P = gaussian(data(labels==i,:), mu(i, :), S);

        L = L + log(sum(pi(i)*P));
    end
    
end

%L = log(l);

end

