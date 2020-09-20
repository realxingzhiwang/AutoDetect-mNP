function [opt_results, step_results] = iterativeclustering(data, n)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


clustering = P_clustering(data);
step_results.step1 = clustering;
clustering = rmfield(clustering, 'All_P');
sub = true(clustering.opt_K, 1);
new_classes = zeros(size(clustering.classes));
new_sub = sub;
cnt = 2;

while cnt <= n && sum(sub) > 0
    sub_clustering = cell(clustering.opt_K, 1);
    for i = 1:clustering.opt_K
        populations = sum(clustering.classes==i);
        if sub(i) && populations > 10
            sub_clustering{i} = P_clustering(data(clustering.classes==i, :));
            if sub_clustering{i}.loglikelihood> clustering.loglikelihood && sub_clustering{i}.opt_K > 1% && min(populations) > 10 % sub_clustering{i}.P_max < 0.2 && 
                new_sub = [new_sub(1:max(new_classes));...
                    true(sub_clustering{i}.opt_K, 1);...
                    new_sub(max(new_classes)+2:end)];
                new_classes(clustering.classes==i) = ...
                    sub_clustering{i}.classes + max(new_classes);
            else
                new_sub(max(new_classes)+1) = false;
                new_classes(clustering.classes==i) = 1 + max(new_classes);
            end
        else
            new_classes(clustering.classes==i) = 1 + max(new_classes);
        end
    end
    sub = new_sub;
    clustering.classes = new_classes;
    clustering.opt_K = max(new_classes);
    [mu, sigma, pi] = compute_distribution(clustering.data, clustering.classes);
    clustering.loglikelihood = loglikelihood(clustering.data, clustering.classes, mu, sigma, pi);
    P_mat = P_matrix(data, new_classes);
    if length(P_mat)~=1
        clustering.P_max = max(P_mat(~eye(clustering.opt_K)));
    else
        clustering.P_max = 1;
    end
    clustering.sub_classes = sub_clustering;
    step_results.(['step' num2str(cnt)]) = clustering;
    cnt = cnt + 1;
    new_classes = zeros(size(clustering.classes));
end

opt_results = rmfield(clustering, 'sub_classes');


end

