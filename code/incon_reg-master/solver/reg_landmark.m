function source_vertex_reg = reg_landmark(source_face, source_vertex, source_vertex_reg_pre,...
                                landmark_source_index, landmark_target_pos,...
                                source_boundary_index, F2Vm, V2Fm, L, param)

K1 = param.UpperBound; % upper bound
K2 = param.LowerBound; % lower bound
ALPHA = param.alpha; % smoothing alpha
BETA = param.beta; %smoothing step
SMOOTH_ITER = param.smooth_iter;

source_vertex_reg = bdsv(source_vertex_reg_pre, source_vertex, source_face, landmark_source_index,...
              landmark_target_pos, K1, K2);

constraint_index = [source_boundary_index; landmark_source_index'];
boundary_pos = source_vertex_reg(source_boundary_index,1:2);
constraint_pos = [boundary_pos; landmark_target_pos];


% compute solution of the Beltrami equation
nu = compute_bc(source_face, source_vertex, source_vertex_reg, 2);
nu(abs(nu)>1) = 0;
nu = heat_flow(nu, F2Vm, V2Fm, L, SMOOTH_ITER, ALPHA, BETA);
source_vertex_reg = lsqc_solver(source_face, source_vertex, nu, constraint_index, constraint_pos);
                         
                            
end 