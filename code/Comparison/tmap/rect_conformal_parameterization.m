function rect = rect_conformal_parameterization(v,f,corner)

% Compute the rectangular conformal parameterization via a novel 
% combination of disk conformal map and disk-to-rectangle conformal map.
%
% Usage:
% map = rect_conformal_parameterization(v,f,corner)
% 
% Input:
% v: nv x 3 vertex coordinates of a simply-connected open triangle mesh
% f: nf x 3 triangulations of a simply-connected open triangle mesh
% (optional) corner: 4 x 1 vertex indices for the four corners of the 
% rectangle, with anti-clockwise orientation
% 4 - 3
% |   |
% 1 - 2
% 
% Output:
% map: nv x 2 vertex coordinates of the rectangular conformal parameterization
% 
% Remarks:
% 1. Please make sure that the input mesh does not contain any 
%    unreferenced vertices/non-manifold vertices/non-manifold edges.
% 2. Please remove all valence 1 boundary vertices (i.e. vertices with 
%    only 1 face attached to them) before running the program.
% 3. Please make sure that the input triangulations f are with 
%    anti-clockwise orientation.
% 4. The output rectangular domain will always have width = 1, while the
%    height depends on the choice of the corners and may not be 1.
%    (The Riemann mapping theorem guarantees that there exists a conformal  
%    map from any simple-connected open surface to the unit square, but if 
%    four vertices on the surface boundary are specified to be the four 
%    corners of the planar domain, the theorem is no longer applicable.)
% 
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi, H. L. Chan, R. Yong, S. Ranjitkar, A. Brook, G. Townsend, K. Chen, and L. M. Lui, 
%     "Tooth morphometry using quasi-conformal theory."
%     Pattern Recognition, 99, 107064, 2020.
%
% Copyright (c) 2019, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

%% Step 1: Linear disk conformal parameterization (Choi, Lui, ACOM 2018)

disk = lineardiskmap(v,f); % Linear disk map (Choi, Lui, ACOM 2018)

%% Step 2: Disk-to-rectangle conformal map (Meng, Choi, Lui, SIIMS 2016)

% Map it to square first
rect = Polygon2Square(disk,f,corner);
sol_u = rect(:,1);
sol_v = rect(:,2);

% Adjust the height of the square to achieve conformal map
mu = 0.*f(:,1);
E_handle = @(h) h_energy(h,sol_u,sol_v,v,f,mu);
h_opt = fminbnd(E_handle,0,10);
rect = [sol_u, h_opt*sol_v];

end

function E = h_energy( h,sol_u,sol_v,v,f,mu)
h_mu = beltrami_coefficient([sol_u,h*sol_v],f,v);
E = sqrt(sum(abs(h_mu(:,1)-mu).^2));
end
