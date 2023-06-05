function [output_struct] = kstest_for_classification_results(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin==1
    data_struct1 = varargin{1};
    data_struct2 = [];
elseif nargin==2
    data_struct1 = varargin{1};
    data_struct2 = varargin{2};
end

output_struct = struct;

dist1 = cell(1, max(data_struct1.classification));
for i=1:length(dist1)
    dist1{i} = data_struct1.features(data_struct1.classification==i);
end

p1 = zeros(length(dist1));
for n=1:size(p1, 1)
    for m=1:size(p1, 2)
        [~, p1(n, m)] = kstest2(dist1{n}, dist1{m});
    end
end

output_struct.P1 = p1;

if ~isempty(data_struct2)

    dist2 = cell(2, max(data_struct2.classification));
    for i=1:length(dist2)
        dist2{i} = data_struct2.features(data_struct2.classification==i);
    end
    
    p2 = zeros(length(dist2));
    for n=1:size(p2, 1)
        for m=1:size(p2, 2)
            [~, p2(n, m)] = kstest2(dist2{n}, dist2{m});
        end
    end

    output_struct.P2 = p2;

    p_co = zeros(length(dist1), length(dist2));
    for j=1:size(p_co, 1)
        for k=1:size(p_co, 2)
            [~, p_co(j,k)] = kstest2(dist1{j}, dist2{k});
        end
    end
    
    output_struct.P_co = p_co;

end

end