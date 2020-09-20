function [varargout] = naivebayes(data, initial_labels)
%[labels, params] = naivebayes(data, initial_labels)
%   naivebayes uses a naive Bayes classifier to perform classification of
%   data specified in the input. data is a N-by-p matrix, where N is the
%   number of data points and p is the number of features. Features should
%   be normalized before classification. initial_labels is a N-by-1 vector,
%   with the n-th element being the class label assigned to the n-th
%   data point. labels is an N-by-1 vector of labels assigned to each data
%   point by the naive Bayes classifier. params is a structure array with
%   three fields, mu, sigma, and pi. params.mu is a k-by-p vector, where k
%   is the number of classes, with n-th row being the means of the p
%   features of the n-th datapoint. params.sigma is a k-by-p vector of the
%   standard deviations of each data point. params.pi is a k-by-1 vector
%   with params.pi(i) = P(Ci), where Ci is the i-th class. sum(params.pi) =
%   1.

[mu, sigma, pi] = compute_distribution(data, initial_labels);

L = loglikelihood(data, initial_labels, mu, sigma, pi);

L_new = 0;
count = 1;

while count<1e3
    
    labels = assignlabels(data, mu, sigma);
    [mu, sigma, pi] = compute_distribution(data, labels);
    L_new = loglikelihood(data, labels, mu, sigma, pi);
    if norm(L_new-L)<1e-5 & accumarray(labels, 1)>1 %Stops when loglikelihood L converges. Does not allow classes containing only 1 data point
        break
    else
        L = L_new;
    end
    count = count+1;

end




if count>=1e3
    disp('Maximum number of iterations reached.')
end

params = struct;
params.mu = mu;
params.sigma = sigma;
params.pi = pi;

if nargout ~= 0
    varargout(1) = {labels};
end

if nargout > 1
    varargout(2) = {params};
end


end

