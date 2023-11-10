% Run a pairwise mapping between two molars using the inconsistent shape
% registration method and compare with Procrustes, TPS, ICP, non-rigid ICP
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
data = dir('HDM 2016_preprocessed/*.mat');
data_list = natsortfiles({data.name})';
prefix = 'HDM 2016_preprocessed/';

param.UpperBound = 1.5;
param.LowerBound = 0.8;
param.alpha = 0.01;
param.beta = 0.1;
param.smooth_iter = 3;
param.intensity_iter = 1;
param.demons_iter = 1; % turn higher if you want stronger intensity matching
param.demons_stepsize = 5; % turn higher if you want stronger intensity matching
landmark_iter = 1;
overall_iter = 20;

for i = 5
    %% load source data and preprocess
    disp(['..................i = ',num2str(i)]);
    
    source_mat =  load([prefix,data_list{i}]);
    source_face = source_mat.f;
    source_vertex = source_mat.v;
    source_gauss_curvature = source_mat.K;
    source_mean_curvature = source_mat.H;
    landmark_source_index = source_mat.lm;
    
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
%         disp('furtherst two points in the target will be used for conformal flattening \n');
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
%         disp('Conformal flattening done\n');
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

%% plot the conformal flattening results
figure; patch('Faces',source_face,'Vertices',flat_source_vertex,'FaceColor','w',...
        'EdgeColor',[255,118,117]/255,'LineWidth',0.5);
axis equal tight
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
hold on;
plot(flat_source_vertex(landmark_source_index,1),flat_source_vertex(landmark_source_index,2), '.', 'Color', [214 48 49]/255, 'MarkerSize',30);
axis tight equal off


figure; patch('Faces',target_face,'Vertices',flat_target_vertex,'FaceColor','w',...
        'EdgeColor',[116,185,255]/255,'LineWidth',0.5);
axis equal tight
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
hold on;
hold on;
plot(flat_target_vertex(landmark_target_index,1),flat_target_vertex(landmark_target_index,2),  '.', 'Color', [9 132 227]/255, 'MarkerSize',30);
axis tight equal off

%% plot the planar mapping result
figure;hold on;
gpp_plot_mesh(source_face, source_vertex_reg, 'FaceColor', 'None', 'EdgeColor', [255,118,117]/255,'LineWidth',0.5); 
gpp_plot_mesh(target_face, flat_target_vertex, 'FaceColor', 'None', 'EdgeColor', [116,185,255]/255,'LineWidth',0.5);

plot(source_vertex_reg(landmark_source_index,1),source_vertex_reg(landmark_source_index,2), '.', 'Color', [214 48 49]/255, 'MarkerSize',30);
axis tight equal off

%% compute the 3D registration result
F1 = TriScatteredInterp(flat_target_vertex(:,1),flat_target_vertex(:,2),target_vertex(:,1));
F2 = TriScatteredInterp(flat_target_vertex(:,1),flat_target_vertex(:,2),target_vertex(:,2));
F3 = TriScatteredInterp(flat_target_vertex(:,1),flat_target_vertex(:,2),target_vertex(:,3));

map = zeros(length(source_vertex_reg),3);
map(:,1) = F1(source_vertex_reg(:,1),source_vertex_reg(:,2));
map(:,2) = F2(source_vertex_reg(:,1),source_vertex_reg(:,2));
map(:,3) = F3(source_vertex_reg(:,1),source_vertex_reg(:,2));

%% plot the 3D registration result
figure;hold on;
gpp_plot_mesh(source_face_reg, map, 'FaceColor', 'w', 'EdgeColor', [255,118,117]/255,'LineWidth',0.5); 
gpp_plot_mesh(target_face, target_vertex, 'FaceColor', 'None', 'EdgeColor', [116,185,255]/255,'LineWidth',0.5);
axis equal tight
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
hold on;
plot3(map(landmark_source_index,1),map(landmark_source_index,2),map(landmark_source_index,3),  '.', 'Color', [214 48 49]/255, 'MarkerSize',30);
axis tight equal off
view([-24 28])

%% plot the normalized curvature (mean, gauss)
cmap = flipud(lbmap(length(source_vertex),'RedBlue'));
% figure; patch('Faces',source_face,'Vertices',source_vertex,'FaceColor','interp','FaceVertexCData',source_mean_curvature_normalized, 'EdgeColor','None');
figure; patch('Faces',source_face,'Vertices',source_vertex,'FaceColor','interp','FaceVertexCData',source_gauss_curvature_normalized, 'EdgeColor','None');
colormap(cmap);
colorbar
axis equal tight off
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
view([-24 28])
caxis([-1 1])
set(gca, 'FontSize',18);

cmap = flipud(lbmap(length(target_vertex),'RedBlue'));
% figure; patch('Faces',target_face,'Vertices',target_vertex,'FaceColor','interp','FaceVertexCData',target_mean_curvature_normalized, 'EdgeColor','None');
figure; patch('Faces',target_face,'Vertices',target_vertex,'FaceColor','interp','FaceVertexCData',target_gauss_curvature_normalized, 'EdgeColor','None');
colormap(cmap);
colorbar
axis equal tight off
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
view([-24 28])
caxis([-1 1])
set(gca, 'FontSize',18);

%% plot BC and curvature difference on the registration result
cmap = flipud(lbmap(length(map),'RedBlue'));
figure; patch('Faces',source_face_reg,'Vertices',map,'FaceColor','interp','FaceVertexCData',(F2Vm*(abs(compute_bc(source_face, flat_source_vertex, source_vertex_reg, 2)))), 'EdgeColor','None');
% figure; patch('Faces',source_face_reg,'Vertices',map,'FaceColor','interp','FaceVertexCData',abs(source_mean_curvature_normalized - target_mean_curvature_normalized_reg), 'EdgeColor','None');
% figure; patch('Faces',source_face_reg,'Vertices',map,'FaceColor','interp','FaceVertexCData',abs(source_gauss_curvature_normalized - target_gauss_curvature_normalized_reg), 'EdgeColor','None');

colormap(cmap);
colorbar
axis equal tight off
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
view([-24 28])
set(gca, 'FontSize',18);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison: Procrustes superimposition
[~,~,transform] = procrustes(target_vertex(landmark_target_index,:),source_vertex(landmark_source_index,:));
c = transform.c;
T = transform.T;
b = transform.b;
map_procrustes = b*source_vertex*T; 
map_procrustes(:,1) = map_procrustes(:,1)+c(1,1); 
map_procrustes(:,2) = map_procrustes(:,2)+c(1,2); 
map_procrustes(:,3) = map_procrustes(:,3)+c(1,3); 

figure;hold on;
gpp_plot_mesh(source_face, map_procrustes, 'FaceColor', 'w', 'EdgeColor', [255,118,117]/255,'LineWidth',0.5); 
gpp_plot_mesh(target_face, target_vertex, 'FaceColor', 'w', 'EdgeColor', [116,185,255]/255,'LineWidth',0.5);
axis equal tight
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
hold on;
plot3(map_procrustes(landmark_source_index,1),map_procrustes(landmark_source_index,2),map_procrustes(landmark_source_index,3),  '.', 'Color', [214 48 49]/255, 'MarkerSize',30);
plot3(target_vertex(landmark_target_index,1),target_vertex(landmark_target_index,2),target_vertex(landmark_target_index,3),  '.', 'Color', [9 132 227]/255, 'MarkerSize',30);
axis tight equal off
view([-24 28])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison: ICP
addpath('code/Comparison')
[TR, TT, ER] = icp_matlabcentral(target_vertex',source_vertex',10);
map_icp = (TR * source_vertex' + TT)';

figure;hold on;
gpp_plot_mesh(source_face, map_icp, 'FaceColor', 'w', 'EdgeColor', [255,118,117]/255,'LineWidth',0.5); 
gpp_plot_mesh(target_face, target_vertex, 'FaceColor', 'w', 'EdgeColor', [116,185,255]/255,'LineWidth',0.5);
axis equal tight
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
hold on;
plot3(map_icp(landmark_source_index,1),map_icp(landmark_source_index,2),map_icp(landmark_source_index,3),  '.', 'Color', [214 48 49]/255, 'MarkerSize',30);
plot3(target_vertex(landmark_target_index,1),target_vertex(landmark_target_index,2),target_vertex(landmark_target_index,3),  '.', 'Color', [9 132 227]/255, 'MarkerSize',30);
axis tight equal off
view([-24 28])
rmpath('code/Comparison')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison: Thin plate spline
addpath('code/Comparison')
map_tps = TPS3D(source_vertex(landmark_source_index,:),target_vertex(landmark_target_index,:),source_vertex);

figure;hold on;
gpp_plot_mesh(source_face, map_tps, 'FaceColor', 'w', 'EdgeColor', [255,118,117]/255,'LineWidth',0.5); 
gpp_plot_mesh(target_face, target_vertex, 'FaceColor', 'w', 'EdgeColor', [116,185,255]/255,'LineWidth',0.5);
axis equal tight
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
hold on;
plot3(map_tps(landmark_source_index,1),map_tps(landmark_source_index,2),map_tps(landmark_source_index,3),  '.', 'Color', [214 48 49]/255, 'MarkerSize',30);
% plot3(target_vertex(landmark_target_index,1),target_vertex(landmark_target_index,2),target_vertex(landmark_target_index,3),  '.', 'Color', [9 132 227]/255, 'MarkerSize',30);
axis tight equal off
view([-24 28])
rmpath('code/Comparison')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Comparison: Non-rigid ICP 
addpath(genpath('code/Comparison/nonrigidICP'))
addpath(genpath('code/Comparison/geom3d'))
tic;
[map_nricp,~,~]=nonrigidICPv2(target_vertex,source_vertex,target_face,source_face,10,1);
toc;

figure;hold on;
gpp_plot_mesh(source_face, map_nricp, 'FaceColor', 'w', 'EdgeColor', [255,118,117]/255,'LineWidth',0.5); 
gpp_plot_mesh(target_face, target_vertex, 'FaceColor', 'w', 'EdgeColor', [116,185,255]/255,'LineWidth',0.5);
axis equal tight
ax = gca;               % get the current axis
ax.Clipping = 'off';    % turn clipping off
hold on;
plot3(map_nricp(landmark_source_index,1),map_nricp(landmark_source_index,2),map_nricp(landmark_source_index,3),  '.', 'Color', [214 48 49]/255, 'MarkerSize',30);
plot3(target_vertex(landmark_target_index,1),target_vertex(landmark_target_index,2),target_vertex(landmark_target_index,3),  '.', 'Color', [9 132 227]/255, 'MarkerSize',30);
axis tight equal off
view([-24 28])

rmpath(genpath('code/Comparison/nonrigidICP'))
rmpath(genpath('code/Comparison/geom3d'))