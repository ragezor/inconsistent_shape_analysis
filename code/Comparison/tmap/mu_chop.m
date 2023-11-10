function Smooth_mu = mu_chop(mu,bound,target)
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi, H. L. Chan, R. Yong, S. Ranjitkar, A. Brook, G. Townsend, K. Chen, and L. M. Lui, 
%     "Tooth morphometry using quasi-conformal theory."
%     Pattern Recognition, 99, 107064, 2020.
%
% Copyright (c) 2019, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi
Smooth_mu = mu;
for i = 1:size(mu,1)
    if(abs(mu(i,1)) >= bound )
        Smooth_mu(i,1) = target*(cos(angle(mu(i,1)))+sqrt(-1)*sin(angle(mu(i,1))));
    end
end