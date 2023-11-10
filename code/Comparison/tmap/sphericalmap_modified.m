function vmap = sphericalmap_modified(v,f,v_ori_new,M)

% If you use this code in your own work, please cite the following paper:
% [1] G. P.-T. Choi and L. M. Lui, 
%     "A Linear Formulation for Disk Conformal Parameterization of Simply-Connected Open Surfaces."
%     Advances in Computational Mathematics, 44(1), pp. 87-114, 2018.
% 
% Copyright (c) 2015-2019, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

%% Find most regular triangle
temp = v(reshape(f',1,length(f)*3),1:3);
e1 = sqrt(sum((temp(2:3:end,1:3) - temp(3:3:end,1:3))'.^2))';
e2 = sqrt(sum((temp(1:3:end,1:3) - temp(3:3:end,1:3))'.^2))';
e3 = sqrt(sum((temp(1:3:end,1:3) - temp(2:3:end,1:3))'.^2))';
regularity = abs(e1./(e1+e2+e3)-1/3)+...
    abs(e2./(e1+e2+e3)-1/3)+abs(e3./(e1+e2+e3)-1/3);

[~,bigtri] = min(regularity);

%% Map the mesh to big triangle

numofv = size(v,1); 
p1 = f(bigtri,1);
p2 = f(bigtri,2);
p3 = f(bigtri,3);

fixed = [p1,p2,p3];
[mrow,mcol,mval] = find(M(fixed,:));
M = M - sparse(fixed(mrow),mcol,mval,numofv, numofv) + ...
        sparse(fixed,fixed,[1,1,1],numofv, numofv);

% Set boundary condition for big triangle
x1 = 0; y1 = 0; x2 = 1; y2 = 0;
a = v(p2,1:3) - v(p1,1:3);
b = v(p3,1:3) - v(p1,1:3);
sin1 = (norm(cross(a,b),2))/(norm(a,2)*norm(b,2));
ori_h = norm(b,2)*sin1;
ratio = norm([x1-x2,y1-y2],2)/norm(a,2);
y3 = ori_h*ratio;
x3 = sqrt(norm(b,2)^2*ratio^2-y3^2);

% Solve matrix equation
c = zeros(numofv,1); c(p1) = x1; c(p2) = x2; c(p3) = x3;
d = zeros(numofv,1); d(p1) = y1; d(p2) = y2; d(p3) = y3;
z = M \ complex(c,d);
z = z-mean(z);
S = [2*real(z)./(1+abs(z).^2), 2*imag(z)./(1+abs(z).^2), (-1+abs(z).^2)./(1+abs(z).^2)];

% Find optimal big triangle size
w = complex(S(:,1)./(1+S(:,3)), S(:,2)./(1+S(:,3)));
[~, index] = sort(abs(z(f(:,1)))+abs(z(f(:,2)))+abs(z(f(:,3))));
inner = index(1);
if inner == bigtri
    inner = index(2);
end

NorthTriSide = (abs(z(f(bigtri,1))-z(f(bigtri,2))) + ...
    abs(z(f(bigtri,2))-z(f(bigtri,3))) + ...
    abs(z(f(bigtri,3))-z(f(bigtri,1))))/3;
SouthTriSide = (abs(w(f(inner,1))-w(f(inner,2))) + ...
    abs(w(f(inner,2))-w(f(inner,3))) + ...
    abs(w(f(inner,3))-w(f(inner,1))))/3;

z = z*(sqrt(NorthTriSide*SouthTriSide))/(NorthTriSide);
S = [2*real(z)./(1+abs(z).^2), -2*imag(z)./(1+abs(z).^2), (-1+abs(z).^2)./(1+abs(z).^2)];

if sum(sum(isnan(S))) ~= 0
    % use tutte map as the initial map
    S = spherical_tutte_map(f,bigtri);
end

%% 
[~,I] = sort(S(:,3));
% number of fixed points needed depends on mesh size, just set it to be 50
fixed = I(1:min(length(v),50)); 
fixnum = length(fixed);

P = [S(:,1)./(1+S(:,3)), S(:,2)./(1+S(:,3))];

%% South pole LBS
mu = beltrami_coefficient(P, f, v_ori_new);
map = linear_beltrami_solver(P,f,mu,fixed,P(fixed,1:2));

% if the result has NaN entries, then most probably the number of boundary constraints is not large enough
while sum(sum(isnan(map))) ~= 0
    % increase the number of boundary constrains and run again
    fixnum = fixnum*5;
    fixed = I(1:min(length(v),fixnum)); 
    map = linear_beltrami_solver(P,f,mu,fixed,P(fixed,1:2)); 
end

z = complex(map(:,1),map(:,2));

vmap = [2*real(z)./(1+abs(z).^2), 2*imag(z)./(1+abs(z).^2), -(abs(z).^2-1)./(1+abs(z).^2)];
