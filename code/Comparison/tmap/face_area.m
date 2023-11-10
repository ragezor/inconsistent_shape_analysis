function fa = face_area(vertex,face)
% If you use this code in your own work, please cite the following paper:
% [1] G. P. T. Choi, H. L. Chan, R. Yong, S. Ranjitkar, A. Brook, G. Townsend, K. Chen, and L. M. Lui, 
%     "Tooth morphometry using quasi-conformal theory."
%     Pattern Recognition, 99, 107064, 2020.
%
% Copyright (c) 2019, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

fi = face(:,1);
fj = face(:,2);
fk = face(:,3);
vij = vertex(fj,:)-vertex(fi,:);
vjk = vertex(fk,:)-vertex(fj,:);
vki = vertex(fi,:)-vertex(fk,:);
a = sqrt(dot(vij,vij,2));
b = sqrt(dot(vjk,vjk,2));
c = sqrt(dot(vki,vki,2));
s = (a+b+c)/2.0;
fa = sqrt(s.*(s-a).*(s-b).*(s-c));
