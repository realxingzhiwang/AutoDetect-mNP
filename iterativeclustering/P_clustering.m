function [results] = P_clustering(data)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

N_ite = 10;
K_range = 1:10;

num_class = zeros(length(K_range), N_ite);
P_mat = cell(length(K_range), N_ite);
P_max = ones(length(K_range), N_ite);
class_indices = cell(length(K_range), N_ite);

parfor k = K_range
    for cnt = 1:N_ite
        class_idx_temp = kmeans(data, k);
        [class_idx_em_temp, params] =...
            naivebayes(data, class_idx_temp);
        num_class(k, cnt) = length(unique(class_idx_em_temp));
        P_mat{k, cnt} = P_matrix(data, class_idx_em_temp, params);
%         P_mat{k, cnt} = zeros(max(class_idx_em_temp));
%         for i = 1:max(class_idx_em_temp)
%             Pii = mean(mvnpdf(data(class_idx_em_temp==i,:),...
%                 params.mu(i, :), diag(params.sigma(i, :))));
%             for j = 1:max(class_idx_em_temp)
%                 Pij = mean(mvnpdf(data(class_idx_em_temp==i,:),...
%                     params.mu(j, :), diag(params.sigma(j, :))));
%                 P_mat{k, cnt}(i, j) = Pij/Pii;
%             end
%         end
% 
        if max(class_idx_em_temp)~=1
            P_max(k, cnt) = max(P_mat{k, cnt}(~eye(max(class_idx_em_temp))));
        else
            P_max(k, cnt) = 1;
        end
        class_indices{k, cnt} = class_idx_em_temp;

    end
end

results = struct;
results.data = data;
results.All_P = min(P_max, [], 2);
[results.P_max, results.opt_K] = min(results.All_P);
[~, opt_run] = min(P_max(results.opt_K, :));
results.classes = class_indices{results.opt_K, opt_run};

[mu, sigma, pi] = compute_distribution(data, results.classes);
results.loglikelihood = loglikelihood(data, results.classes, mu, sigma, pi);


end

