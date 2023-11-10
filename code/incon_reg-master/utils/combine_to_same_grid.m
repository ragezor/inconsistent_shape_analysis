function [g3,g4,Mg,Sg,mask,scale] = combine_to_same_grid(source_vertex, source_intensity, ...
                                target_vertex, target_intensity,...
                                source_boundary_index, target_boundary_index)
%output Mg and Sg defined one the grid [g3,g4]
%ready to be used for itensity matching

% the rectangular grid (g3, g4) should contain V_t
sx = max(target_vertex(:,1)) - min(target_vertex(:,1));
sy = max(target_vertex(:,2)) - min(target_vertex(:,2));
r = sx/sy; N = round(sqrt(size(target_vertex,1)));
[g4,g3] = ndgrid(1*N:-0.4:0, 0:0.4:1*round(r*N)); %natural ordering

%scale and translate the grid
g4 = g4/N*sy; g3 = g3/r/N*sx;
g4 = g4 + min(target_vertex(:,2)); 
g3 = g3 + min(target_vertex(:,1));

scale = abs((g4(end) - g4(1)) / size(g4,2));
Ftmp = delaunay(g3(:),g4(:));
L = cotmatrix([g3(:),g4(:)],Ftmp);
% interpolating to [g3,g4]
intpM = scatteredInterpolant(source_vertex(:,1),source_vertex(:,2),source_intensity);
intpS= scatteredInterpolant(target_vertex(:,1),target_vertex(:,2),target_intensity);
Mg = intpM(g3,g4);
Sg = intpS(g3,g4);
bdy_xym = source_vertex(source_boundary_index,:);
bdy_xyt = target_vertex(target_boundary_index,:);
bout_m = ~in_or_out(bdy_xym,[g3(:),g4(:)]); bout_m = reshape(bout_m,size(g3,1),size(g3,2));
bout_t = ~in_or_out(bdy_xyt,[g3(:),g4(:)]); bout_t = reshape(bout_t,size(g3,1),size(g3,2));
mask = isnan(Mg) | isnan(Sg) | bout_m | bout_t;
Mg(mask) = 0; Sg(mask) = 0;
mask = ~mask(:)*1.0;
for kkk = 1:10
    mask = mask + 0.1*L*mask;
end
mask = reshape(mask,size(g3,1),size(g3,2)); 
% figure(10); imshow(mask); title('Overlap domain');
Mg = Mg.*mask; Sg = Sg.*mask;
end