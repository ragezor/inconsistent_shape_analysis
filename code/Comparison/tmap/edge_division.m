function [sort_edge,sort_point] = edge_division(B,point)
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi, H. L. Chan, R. Yong, S. Ranjitkar, A. Brook, G. Townsend, K. Chen, and L. M. Lui, 
%     "Tooth morphometry using quasi-conformal theory."
%     Pattern Recognition, 99, 107064, 2020.
%
% Copyright (c) 2019, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi
location = [];
for k = 1:length(point)
    for i = 1:length(B)
        if B(i,1) == point(k)
            location = [location;i];
        end
    end
end

[sort_location,index] = sort(location);
sort_point = point(index);
for i = 1:length(sort_location)-1
    edge{i} = B(sort_location(i,1):sort_location(i+1,1),1);
end
edge{length(sort_location)} = [B(sort_location(end,1):end,1);B(1:sort_location(1,1),1)];

for i = 1:length(point)-1
    for j = 1:length(edge)
        temp = edge{j};
        if sum(ismember([point(i);point(i+1)],[temp(1);temp(end)])) == 2
            sort_edge{i} = temp;
        end
    end
end
for j = 1:length(edge)
    temp = edge{j};
    if sum(ismember([point(end);point(1)],[temp(1);temp(end)])) == 2
        sort_edge{length(point)} = temp;
    end
end
end