function [source_vertex_reg_intersect_index, correspondence_mask, ...
    source_vertex_reg_3D, source_face_reg, displace, dist,...
    target_intensity_reg, intensity_diff] = ...
    prepare_result(source_face, source_vertex, source_vertex_reg,...
                    target_face, target_vertex, flat_target_vertex,...
                    source_intensity, target_intensity)
               
%% find the common domain submesh from inconsistent shape registration

% assume only one boundary component
[~, target_boudary_index] = meshboundaries(target_face);
target_boundary_pos = flat_target_vertex(target_boudary_index,:);

% find the vertices in the common domain
source_vertex_reg_intersect_index = in_or_out(target_boundary_pos, source_vertex_reg);   
source_vertex_reg_intersect_pos = source_vertex_reg(source_vertex_reg_intersect_index,:);
correspondence_mask = zeros(size(source_vertex_reg,1), 1);
correspondence_mask(source_vertex_reg_intersect_index, 1) = 1;
tmp = 1:size(source_vertex, 1);
source_vertex_reg_intersect_index = tmp(source_vertex_reg_intersect_index);
% get barycentric coordinate of registration with respect to the target mesh 
Tri_target = triangulation(target_face, flat_target_vertex);
tid = pointLocation(Tri_target, source_vertex_reg_intersect_pos);
source_vertex_reg_bary = cartesianToBarycentric(Tri_target, tid, source_vertex_reg_intersect_pos);
% use the barycentric coordinate to get 3D registration
Tri_target_3D = triangulation(target_face, target_vertex);
source_vertex_reg_3D_sub = barycentricToCartesian(Tri_target_3D, tid, source_vertex_reg_bary);
source_vertex_reg_3D = arap(source_vertex, source_face, source_vertex_reg_intersect_index', source_vertex_reg_3D_sub);
% find the faces in the common domain
f_list = ismember(source_face(:,1),source_vertex_reg_intersect_index') | ismember(source_face(:,2),source_vertex_reg_intersect_index') |...
    ismember(source_face(:,3), source_vertex_reg_intersect_index');
source_face_reg = source_face(f_list,:);
% interpolate target intensity to registered
target_intensity_reg = griddata(flat_target_vertex(:,1), flat_target_vertex(:,2), target_intensity, source_vertex_reg(:,1), source_vertex_reg(:,2));
target_intensity_reg(isnan(target_intensity_reg)) = 0;
% distance from initial alignment
displace = source_vertex_reg_3D - source_vertex;
% displace_3D = displace_3D - mean(displace_3D, 2);
dist = sqrt(sum((displace).^2, 2));
intensity_diff = abs(source_intensity - target_intensity_reg);
end