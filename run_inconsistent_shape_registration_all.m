% run inconisistent shape registration for all 50x50 tooth meshes
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
prefix = 'HDM 2016_preprocessed/';
data = dir([prefix,'*.mat']);
data_list = natsortfiles({data.name})';

% parameters
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

Hndiff_all = zeros(length(data_list),length(data_list)); 
Kndiff_all = zeros(length(data_list),length(data_list)); 
mu_all = zeros(length(data_list),length(data_list));

%% all 50x50 mappings, take several hours
for i = 1:length(data_list)
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
    parfor j = 1:length(data_list)
        try
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
            
            distance_table = repmat(reshape(target_vertex(landmark_target_index, :), ...
                [1, num_landmark, 3]), [num_landmark, 1, 1]) - ...
                repmat(reshape(target_vertex(landmark_target_index, :), ...
                [num_landmark, 1, 3]), [1, num_landmark, 1]);
            distance_table = sum(distance_table.^2, 3);
            
            [~, sorted_] = sort(distance_table(:));
            lscm_target_ind1 = mod(sorted_(end), num_landmark);
            lscm_target_ind2 = round(sorted_(end)/num_landmark + 0.5);
            lscm_target_ind = [landmark_target_index(lscm_target_ind1), landmark_target_index(lscm_target_ind2)];
            lscm_source_ind = [landmark_source_index(lscm_target_ind1), landmark_source_index(lscm_target_ind2)];
            %% conformal flattening
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

            if initial_intensity_error < 1e-10 && initial_landmark_error < 1e-10
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
            Hndiff_all(i,j) = sum(source_vertex_area(source_vertex_reg_intersect_index).*Hndiff)/sum(source_vertex_area(source_vertex_reg_intersect_index));  
            Kndiff_all(i,j) = sum(source_vertex_area(source_vertex_reg_intersect_index).*Kndiff)/sum(source_vertex_area(source_vertex_reg_intersect_index));  
            mu_all(i,j) = sum(source_face_area(f_list).*mu)/sum(source_face_area(f_list));  
        
        catch
            
        end
    end
end

save('inconsistent_surface_registration_result.mat','data_list','mu_all','Hndiff_all','Kndiff_all');

