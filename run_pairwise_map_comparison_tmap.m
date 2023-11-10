% Run a pairwise mapping between two inconsistently segmented molars using 
% the inconsistent shape registration method and compare with global T-map
% 
% Reference:
% G. P. T. Choi, D. Qiu and L. M. Lui, 
% "Shape analysis via inconsistent surface registration."
% Proceedings of the Royal Society A, 476(2242), 20200147, 2020.
% 
% Copyright (c) 2020, Gary Pui-Tung Choi, Di Qiu, Lok Ming Lui

addpath('code')
addpath('code/natsortfiles')
addpath(genpath('code/incon_reg-master'))
addpath(genpath('code/GPToolBox-master'))

%%
param.UpperBound = 1.5;
param.LowerBound = 0.8;
param.alpha = 0.01;
param.beta = 0.1;
param.smooth_iter = 3;
param.intensity_iter = 1;
param.demons_iter = 1; 
param.demons_stepsize = 5;
landmark_iter = 1;
overall_iter = 20;

for i = 5
    %% load source data and preprocess
    
    source_mat =  load([prefix,data_list{i}]);
    source_face_ori = source_mat.f;
    source_vertex_ori = source_mat.v;
    source_gauss_curvature_ori = source_mat.K;
    source_mean_curvature_ori = source_mat.H;
    landmark_source_index_ori = source_mat.lm;
    
    
    % segment the source mesh inconsistently
    id = find((source_vertex_ori(source_face_ori(:,1),3)+source_vertex_ori(source_face_ori(:,2),3)+source_vertex_ori(source_face_ori(:,3),3))/3 <-0.2*10^(-3));
    
    % remove all degree-2 vertices
    f = source_face_ori(setdiff(1:length(source_face_ori), id),:);
    A = sparse([f(:,1);f(:,2);f(:,3);f(:,1);f(:,2);f(:,3)],[f(:,2);f(:,3);f(:,1);f(:,3);f(:,1);f(:,2)],ones(length(f)*6,1))/2;
    A = ceil(A);
    vid_degree_2 = find(sum(A,2) == 2);
    fid_degree_2 = zeros(size(vid_degree_2));
    for iii = 1:length(vid_degree_2)
        fid_degree_2_temp = find(f(:,1) == vid_degree_2(iii) | f(:,2) == vid_degree_2(iii) | f(:,3) == vid_degree_2(iii));
        fid_degree_2(iii) = find(source_face_ori(:,1) == f(fid_degree_2_temp,1) & source_face_ori(:,2) == f(fid_degree_2_temp,2) & source_face_ori(:,3) == f(fid_degree_2_temp,3));
    end
    
    f = source_face_ori(setdiff(1:length(source_face_ori), [id;fid_degree_2]),:);
    A = sparse([f(:,1);f(:,2);f(:,3);f(:,1);f(:,2);f(:,3)],[f(:,2);f(:,3);f(:,1);f(:,3);f(:,1);f(:,2)],ones(length(f)*6,1))/2;
    A = ceil(A);
    vid_degree_2 = find(sum(A,2) == 2);
    fid_degree_2_2 = zeros(size(vid_degree_2));
    for iii = 1:length(vid_degree_2)
        fid_degree_2_temp = find(f(:,1) == vid_degree_2(iii) | f(:,2) == vid_degree_2(iii) | f(:,3) == vid_degree_2(iii));
        fid_degree_2_2(iii) = find(source_face_ori(:,1) == f(fid_degree_2_temp,1) & source_face_ori(:,2) == f(fid_degree_2_temp,2) & source_face_ori(:,3) == f(fid_degree_2_temp,3));
    end
    
    [source_face, source_vertex, father] = gpp_clean_mesh(source_face_ori(setdiff(1:length(source_face_ori), [id;fid_degree_2; fid_degree_2_2]),:), source_vertex_ori);
    
    source_gauss_curvature = source_gauss_curvature_ori(father);
    source_mean_curvature = source_mean_curvature_ori(father);
    landmark_source_index = landmark_source_index_ori; 
    for iii = 1:length(landmark_source_index_ori)
        landmark_source_index(iii) = find(father == landmark_source_index_ori(iii));
    end
    
    % range: -1 to 1
    source_mean_curvature_normalized = source_mean_curvature/max(abs(source_mean_curvature));
    source_gauss_curvature_normalized = source_gauss_curvature/max(abs(source_gauss_curvature));
    
    num_landmark = length(landmark_source_index);
    
    % build laplacian of the source domain for coefficient smoothing
    L = cotmatrix(source_vertex, source_face);
    F2Vm = F2V(source_vertex', source_face');
    V2Fm = V2F(source_vertex', source_face');
        
    % assume one boundary component                           
    [~, source_boundary_index] = meshboundaries(source_face);
    
    % for evaluating the shape difference
    source_vertex_area = vertex_area(source_face,source_vertex);
    source_face_area = face_area(source_face,source_vertex);
    
    % range: 0 to 1
    source_intensity = intensity_normalization(source_gauss_curvature);
        
    %% Target
    for j = 29
        %% load target data and preprocess
        target_mat = load([prefix,data_list{j}]);

        target_face = target_mat.f;
        target_vertex = target_mat.v;
        target_gauss_curvature = target_mat.K;
        target_mean_curvature = target_mat.H;
        landmark_target_index = target_mat.lm;

        % choose matching intensity to be gauss curvature
        % range: 0 to 1
        target_intensity = intensity_normalization(target_gauss_curvature);
        
        % range: -1 to 1
        target_mean_curvature_normalized = target_mean_curvature/max(abs(target_mean_curvature));
        target_gauss_curvature_normalized = target_gauss_curvature/max(abs(target_gauss_curvature));
        %%
        distance_table = repmat(reshape(target_vertex(landmark_target_index, :), ... 
            [1, num_landmark, 3]), [num_landmark, 1, 1]) - ...
            repmat(reshape(target_vertex(landmark_target_index, :), ...
            [num_landmark, 1, 3]), [1, num_landmark, 1]);
        distance_table = sum(distance_table.^2, 3);
        %%
        [~, sorted_] = sort(distance_table(:));
        lscm_target_ind1 = mod(sorted_(end), num_landmark);
        lscm_target_ind2 = round(sorted_(end)/num_landmark + 0.5);
        lscm_target_ind = [landmark_target_index(lscm_target_ind1), landmark_target_index(lscm_target_ind2)];
        lscm_source_ind = [landmark_source_index(lscm_target_ind1), landmark_source_index(lscm_target_ind2)];
%%
        % conformal flattening
        flat_source_vertex = lscm(source_vertex, source_face, lscm_source_ind, [0, 0; 1, 0]);
        flat_target_vertex = lscm(target_vertex, target_face, lscm_target_ind, [0, 0; 1, 0]);
        landmark_target_pos = flat_target_vertex(landmark_target_index,1:2);
        % assume one boundary component                           
        [~, target_boundary_index] = meshboundaries(target_face);

        %% Inconsistent registration
        initial_landmark_error = compute_landamrk_err(flat_source_vertex, landmark_source_index,...
            landmark_target_pos);
        [~,~,source_intensity_grid, target_intensity_grid, ~,~] = combine_to_same_grid(...
                                                flat_source_vertex, source_intensity,...
                                                flat_target_vertex, target_intensity,...
                                                source_boundary_index, target_boundary_index);
        initial_intensity_error = sum(abs(source_intensity_grid(:) - target_intensity_grid(:))); 

        if initial_intensity_error == 0 && initial_landmark_error == 0
            source_vertex_reg = flat_source_vertex;
        else
            % algo begins
            source_vertex_reg_pre = flat_source_vertex;
            for iter = 1:overall_iter
                fprintf('Iter %d \n', iter);
                [source_vertex_reg, ie] = reg_intensity(source_face, flat_source_vertex, source_vertex_reg_pre,...
                                            source_intensity, flat_target_vertex, target_intensity,...
                                            source_boundary_index, target_boundary_index,...
                                            F2Vm, V2Fm, L, param);

    
                source_vertex_reg = reg_landmark(source_face, flat_source_vertex, source_vertex_reg,...
                                            landmark_source_index, landmark_target_pos,...
                                            source_boundary_index, F2Vm, V2Fm, L, param);

    
                % calculate the landmark and intensity difference
                le = compute_landamrk_err(source_vertex_reg, ...
                                    landmark_source_index,...
                                    landmark_target_pos);
                %
                source_vertex_reg_pre = source_vertex_reg;
                
            end
            % algo ends
        end
        %% Get the results
        [source_vertex_reg_intersect_index, correspondence_mask, ...
            source_vertex_reg_3D, source_face_reg, displace, dist,...
            target_intensity_reg, intensity_diff] = ...
            prepare_result(source_face, source_vertex, source_vertex_reg,...
                            target_face, target_vertex, flat_target_vertex,...
                            source_intensity, target_intensity);
                    
        % interpolate target curvature with normalization to registered
        target_mean_curvature_normalized_reg = griddata(flat_target_vertex(:,1), flat_target_vertex(:,2), target_mean_curvature_normalized, source_vertex_reg(:,1), source_vertex_reg(:,2));
        target_mean_curvature_normalized_reg(isnan(target_mean_curvature_normalized_reg)) = 0;
        target_gauss_curvature_normalized_reg = griddata(flat_target_vertex(:,1), flat_target_vertex(:,2), target_gauss_curvature_normalized, source_vertex_reg(:,1), source_vertex_reg(:,2));
        target_gauss_curvature_normalized_reg(isnan(target_gauss_curvature_normalized_reg)) = 0;
        
        % normalized curvature difference over the common domain
        Hndiff = (abs(source_mean_curvature_normalized(source_vertex_reg_intersect_index) - target_mean_curvature_normalized_reg(source_vertex_reg_intersect_index)));
        Kndiff = (abs(source_gauss_curvature_normalized(source_vertex_reg_intersect_index) - target_gauss_curvature_normalized_reg(source_vertex_reg_intersect_index)));
        
        % BC over the common domain
        f_list = ismember(source_face(:,1),source_vertex_reg_intersect_index') | ismember(source_face(:,2),source_vertex_reg_intersect_index') |...
                    ismember(source_face(:,3), source_vertex_reg_intersect_index');

        mu = abs(compute_bc(source_face, flat_source_vertex, source_vertex_reg, 2));   
        mu = mu(f_list,:);
        
        %% scalar quantities capturing the shape difference
        % integral of those quantities over the overlapping domain, divided by the surface area
        Hndiff = sum(source_vertex_area(source_vertex_reg_intersect_index).*Hndiff)/sum(source_vertex_area(source_vertex_reg_intersect_index));  
        Kndiff = sum(source_vertex_area(source_vertex_reg_intersect_index).*Kndiff)/sum(source_vertex_area(source_vertex_reg_intersect_index));  
        mu = sum(source_face_area(f_list).*mu)/sum(source_face_area(f_list));  
        
    end
end

%% plot source and target mesh
figure; patch('Faces',source_face,'Vertices',source_vertex,'FaceColor','w',...
        'EdgeColor',[255,118,117]/255,'LineWidth',0.5);
axis equal tight
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
hold on;
plot3(source_vertex(landmark_source_index,1),source_vertex(landmark_source_index,2),source_vertex(landmark_source_index,3), '.', 'Color', [214 48 49]/255, 'MarkerSize',30);
axis tight equal off
view([-24 28])

figure; patch('Faces',target_face,'Vertices',target_vertex,'FaceColor','w',...
        'EdgeColor',[116,185,255]/255,'LineWidth',0.5);
axis equal tight
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
hold on;
plot3(target_vertex(landmark_target_index,1),target_vertex(landmark_target_index,2),target_vertex(landmark_target_index,3), '.', 'Color', [9 132 227]/255, 'MarkerSize',30);
axis tight equal off
view([-24 28])


%% compute the 3D registration result
F1 = TriScatteredInterp(flat_target_vertex(:,1),flat_target_vertex(:,2),target_vertex(:,1));
F2 = TriScatteredInterp(flat_target_vertex(:,1),flat_target_vertex(:,2),target_vertex(:,2));
F3 = TriScatteredInterp(flat_target_vertex(:,1),flat_target_vertex(:,2),target_vertex(:,3));

map = zeros(length(source_vertex_reg),3);
map(:,1) = F1(source_vertex_reg(:,1),source_vertex_reg(:,2));
map(:,2) = F2(source_vertex_reg(:,1),source_vertex_reg(:,2));
map(:,3) = F3(source_vertex_reg(:,1),source_vertex_reg(:,2));

% plot the 3D registration result
figure;hold on;
gpp_plot_mesh(target_face, target_vertex, 'FaceColor', 'None', 'EdgeColor', [116,185,255]/255,'LineWidth',0.5);
gpp_plot_mesh(source_face_reg, map, 'FaceColor', 'w', 'EdgeColor', [255,118,117]/255,'LineWidth',0.5); 
axis equal tight
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
hold on;
plot3(map(landmark_source_index,1),map(landmark_source_index,2),map(landmark_source_index,3),  '.', 'Color', [214 48 49]/255, 'MarkerSize',30);
axis tight equal off
view([-24 28])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison: T-map
addpath('code/Comparison/tmap')

v1 = source_vertex;
f1 = source_face;
lm1 = landmark_source_index;

v2 = target_vertex;
f2 = target_face;
lm2 = landmark_target_index;

% Find the four corners of the rectangular domain
bdy1 = meshboundaries(f1); bdy1 = bdy1{1};
corner1 = zeros(4,1);
[~,temp] = min((v1(bdy1,1) - v1(lm1(1),1)).^2 + (v1(bdy1,2) - v1(lm1(1),2)).^2 + (v1(bdy1,3) - v1(lm1(1),3)).^2);
corner1(1) = bdy1(temp);
[~,temp] = min((v1(bdy1,1) - v1(lm1(2),1)).^2 + (v1(bdy1,2) - v1(lm1(2),2)).^2 + (v1(bdy1,3) - v1(lm1(2),3)).^2);
corner1(2) = bdy1(temp);
[~,temp] = min((v1(bdy1,1) - v1(lm1(3),1)).^2 + (v1(bdy1,2) - v1(lm1(3),2)).^2 + (v1(bdy1,3) - v1(lm1(3),3)).^2);
corner1(3) = bdy1(temp);
[~,temp] = min((v1(bdy1,1) - v1(lm1(4),1)).^2 + (v1(bdy1,2) - v1(lm1(4),2)).^2 + (v1(bdy1,3) - v1(lm1(4),3)).^2);
corner1(4) = bdy1(temp);

% Compute the rectangular conformal parameterization
rect1 = rect_conformal_parameterization(v1,f1,corner1);

% TR2 = TriRep(f2,v2); B2 = freeBoundary(TR2); bdy2 = B2(:,1);
bdy2 = meshboundaries(f2); bdy2 = bdy2{1};
corner2 = zeros(4,1);
[~,temp] = min((v2(bdy2,1) - v2(lm2(1),1)).^2 + (v2(bdy2,2) - v2(lm2(1),2)).^2 + (v2(bdy2,3) - v2(lm2(1),3)).^2);
corner2(1) = bdy2(temp);
[~,temp] = min((v2(bdy2,1) - v2(lm2(2),1)).^2 + (v2(bdy2,2) - v2(lm2(2),2)).^2 + (v2(bdy2,3) - v2(lm2(2),3)).^2);
corner2(2) = bdy2(temp);
[~,temp] = min((v2(bdy2,1) - v2(lm2(3),1)).^2 + (v2(bdy2,2) - v2(lm2(3),2)).^2 + (v2(bdy2,3) - v2(lm2(3),3)).^2);
corner2(3) = bdy2(temp);
[~,temp] = min((v2(bdy2,1) - v2(lm2(4),1)).^2 + (v2(bdy2,2) - v2(lm2(4),2)).^2 + (v2(bdy2,3) - v2(lm2(4),3)).^2);
corner2(4) = bdy2(temp);

% Compute the rectangular conformal parameterization
rect2 = rect_conformal_parameterization(v2,f2,corner2);

interior_landmark = landmark_source_index'; % indices of the interior landmarks
interior_landmark_target = rect2(landmark_target_index,1:2); % target positions of the landmarks
height = max(rect2(:,2)); % target height of the rectangle

% Compute the Teichmuller map between the rectangles using QC Iteration
[tmap,tmap_bc] = rectangular_Teichmuller_map(rect1,f1,interior_landmark,interior_landmark_target,corner1,height);

% Computing the Teichmuller distance
d = 1/2*log((1+mean(abs(tmap_bc)))/(1-mean(abs(tmap_bc)))); 
fprintf('Teichmuller distnace between the teeth = %f\n',d);

% Construct the interpolants
Fx = scatteredInterpolant(rect2(:,1),rect2(:,2),v2(:,1),'linear');
Fy = scatteredInterpolant(rect2(:,1),rect2(:,2),v2(:,2),'linear');
Fz = scatteredInterpolant(rect2(:,1),rect2(:,2),v2(:,3),'linear');

% Find the surface T-map
tmap_surface = zeros(length(tmap),3);
tmap_surface(:,1) = Fx(tmap(:,1),tmap(:,2));
tmap_surface(:,2) = Fy(tmap(:,1),tmap(:,2));
tmap_surface(:,3) = Fz(tmap(:,1),tmap(:,2));

% plot the 3D registration result
figure;hold on;
gpp_plot_mesh(target_face, target_vertex, 'FaceColor', 'None', 'EdgeColor', [116,185,255]/255,'LineWidth',0.5);
gpp_plot_mesh(source_face_reg, tmap_surface, 'FaceColor', 'w', 'EdgeColor', [255,118,117]/255,'LineWidth',0.5); 
axis equal tight
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
hold on;
plot3(tmap_surface(landmark_source_index,1),tmap_surface(landmark_source_index,2),tmap_surface(landmark_source_index,3),  '.', 'Color', [214 48 49]/255, 'MarkerSize',30);
axis tight equal off
view([-24 28])

rmpath('code/Comparison/tmap')