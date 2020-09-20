function [P_mat] = P_matrix(data, labels, params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% 
% if isstruct(varargin)
%     data = varargin.data;
%     classification = varargin.classes;
%     params = varargin.params;
% else
%     data = varargin{1};
%     classification = varargin{2};
%     params = varargin{3};
% end


if max(labels) == 1
    P_mat = 1;
    return
end

if nargin < 3
   [params.mu, params.sigma, params.pi] = compute_distribution(data, labels);
end
    
P_mat = zeros(max(labels));
for i = 1:max(labels)
    Pii = mean(mvnpdf(data(labels==i,:),...
        params.mu(i, :), diag(params.sigma(i, :))));
    for j = 1:max(labels)
        Pij = mean(mvnpdf(data(labels==i,:),...
            params.mu(j, :), diag(params.sigma(j, :))));
        P_mat(i, j) = Pij/Pii;
    end
end

end

