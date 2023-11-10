function [source_vertex_reg, intensity_err] = reg_intensity(source_face, source_vertex, source_vertex_reg_pre,...
                                source_intensity, target_vertex, target_intensity,...
                                source_boundary_index, target_boundary_index,...
                                F2Vm, V2Fm, L, param)

ITER = param.intensity_iter;
SUB_ITER = param.demons_iter;
STEP_SIZE = param.demons_stepsize;
SMOOTH_ITER = param.smooth_iter;
ALPHA = param.alpha; % smoothing alpha
BETA = param.beta; %smoothing step

source_vertex_reg = source_vertex_reg_pre;

    
[g1, g2, source_intensity_grid, target_intensity_grid, mask, scale] = combine_to_same_grid(...
                                        source_vertex_reg, source_intensity,...
                                        target_vertex, target_intensity,...
                                        source_boundary_index, target_boundary_index);  
iter = 0;

while iter < ITER
        
        [Txg,Tyg,intensity_err] = intensity_fitting(source_intensity_grid, target_intensity_grid, SUB_ITER, STEP_SIZE, scale);
%         Txg = Txg.*mask; Tyg = Tyg.*mask;
        Tx = griddata(g1(:), g2(:), Txg(:), source_vertex_reg(:,1), source_vertex_reg(:,2));
        Ty = griddata(g1(:), g2(:), Tyg(:), source_vertex_reg(:,1), source_vertex_reg(:,2));
        Tx(isnan(Tx)) = 0; 
        Ty(isnan(Ty)) = 0; 
        boundary_pos = source_vertex_reg + [-Tx,Ty]; 
        boundary_pos = boundary_pos(source_boundary_index, 1:2);
        source_vertex_reg = source_vertex_reg + [-Tx,Ty]; 
                
        nu = compute_bc(source_face, source_vertex, source_vertex_reg, 2);        
        nu(abs(nu)>1) = 0;
%         nu = heat_flow(nu,F2Vm,V2Fm,L,SMOOTH_ITER,ALPHA,BETA);
                
        source_vertex_reg = lsqc_solver(source_face, source_vertex, nu,...
                    source_boundary_index, boundary_pos);  
        iter = iter + 1;
end

% intensity_err = intensity_err/size(source_intensity_grid,1)/size(source_intensity_grid,2);