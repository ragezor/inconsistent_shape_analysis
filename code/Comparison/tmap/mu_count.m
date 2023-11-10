function count = mu_count(mu,bound)
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi, H. L. Chan, R. Yong, S. Ranjitkar, A. Brook, G. Townsend, K. Chen, and L. M. Lui, 
%     "Tooth morphometry using quasi-conformal theory."
%     Pattern Recognition, 99, 107064, 2020.
%
% Copyright (c) 2019, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi
count = 0;
for t = 1:size(mu,1);
    if(abs(mu(t,1)) >= bound)
        count = count + 1;
    end
end
