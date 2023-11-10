function [map,map_mu] = rectangular_Teichmuller_map(v,f,point,target,corner,height)
% Find the landmark-matching Teichmuller map between two rectangular domains
%
% Input:
% v: nv x 2 vertex coordinates of the source rectangular domain
% f: nf x 3 triangulations of the source rectangular domain
% point: k x 1 indices of the landmarks
% target: k x 2 vertex coordinates of the target positions for the landmarks 
% corner: the 4 corners of the rectangle in anticlockwise direction starting from [0,0]
%        i.e., 
%        4 - 3
%        |   |
%        1 - 2
% height: height of the target rectangular domain
% 
% Output:
% map: nv x 2 vertex coordinates of the Teichmuller map
%
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi, H. L. Chan, R. Yong, S. Ranjitkar, A. Brook, G. Townsend, K. Chen, and L. M. Lui, 
%     "Tooth morphometry using quasi-conformal theory."
%     Pattern Recognition, 99, 107064, 2020.
%
% Copyright (c) 2019, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

%%
mu = zeros(length(f),1);
TR = TriRep(f,v(:,1),v(:,2));
ring = vertexAttachments(TR);
B = freeBoundary(TR);
edge = edge_division(B,corner);
temp_edge = edge;
temp_corner = [corner;corner(1)];
for i = 1:4
    for j = 1:4
        if sum(ismember(temp_corner(i:i+1),edge{j})) == 4
            temp_edge{i} = edge{j};
        end
    end
end

vertical = [[edge{4};edge{2}],[zeros(length(edge{4}),1);ones(length(edge{2}),1)]];
horizontal = [[edge{1};edge{3}],[zeros(length(edge{1}),1);height*ones(length(edge{3}),1)]];

if ~isempty(point)
    constraint_x = [vertical;point,target(:,1)];
    constraint_y = [horizontal;point,target(:,2)];
else
    constraint_x = vertical;
    constraint_y = horizontal;
end

setup.iter_tol = .01;
setup.bound = 1;
setup.Chop_val = setup.bound - 0.05;
setup.lm_err = 0;

map = map_reconstruct(mu,v,f,ring,constraint_x,constraint_y);
update_mu = beltrami_coefficient(v,f,map);
mu_Diff = 1;
count = 1;
numofiter = 0;

%% QC Iteration
while (mu_Diff > setup.iter_tol) || (count > 0)
    
    Prev_mu = update_mu;
    update_mu = mu_chop(update_mu,1,0);
    
    smooth_mu = mean(abs(update_mu)).*(cos(angle(real(update_mu)+sqrt(-1)*imag(update_mu)))+...
        sqrt(-1)*sin(angle(real(update_mu)+sqrt(-1)*imag(update_mu))));
    
    map = map_reconstruct(smooth_mu,v,f,ring,constraint_x,constraint_y);
    update_mu = beltrami_coefficient(v,f,map);
    mu_Diff = max(abs(update_mu - Prev_mu));
    count = mu_count(update_mu,setup.bound);
    
    numofiter = numofiter + 1;
    if numofiter > 500
        break;
    end
end
map_mu = beltrami_coefficient(v,f,map);

end
