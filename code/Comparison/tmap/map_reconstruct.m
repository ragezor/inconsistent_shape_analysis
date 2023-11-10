function map = map_reconstruct(mu,v,f,ring,constraint_x,constraint_y)
% Reconstruct the quasi-conformal map with the constraints.
%
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi, H. L. Chan, R. Yong, S. Ranjitkar, A. Brook, G. Townsend, K. Chen, and L. M. Lui, 
%     "Tooth morphometry using quasi-conformal theory."
%     Pattern Recognition, 99, 107064, 2020.
%
% Copyright (c) 2019, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

[af,bf,gf,tri_area] = compute_info(v,f,real(mu),imag(mu));

if length(tri_area(any(tri_area,2),:)) ~= length(tri_area)
    error('Triangles with zero area exists !');
end

[I,J,Val,d1] = defineA(v,f,af,bf,gf,tri_area,constraint_x,ring);
A = sparse(I,J,Val,size(v,1),size(v,1));
umap = A\d1;

[I2,J2,Val2,d2] = defineB(v, f, af,bf,gf,tri_area,constraint_y,ring);
B = sparse(I2,J2,Val2,size(v,1),size(v,1));
vmap = B\d2;

map = [umap,vmap];

end