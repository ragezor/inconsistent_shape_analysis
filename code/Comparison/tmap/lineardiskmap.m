function map = lineardiskmap(v,f)
% Conformally map a simply-connected open triangle mesh to the unit disk using a linear formulation
%
% Input:
% v: nv x 3 vertex coordinates of a simply-connected open triangle mesh
% f: nf x 3 triangulations of a simply-connected open triangle mesh
% 
% Output:
% map: nv x 3 vertex coordinates of the disk conformal map
% 
% Remark:
% 1. Please make sure that the input mesh does not contain any unreferenced vertices/non-manifold vertices/non-manifold edges.
% 2. Please remove all valence 1 boundary vertices (i.e. vertices with only 1 face attached to them) before running the program.
%
% If you use this code in your own work, please cite the following paper:
% [1] G. P.-T. Choi and L. M. Lui, 
%     "A Linear Formulation for Disk Conformal Parameterization of Simply-Connected Open Surfaces."
%     Advances in Computational Mathematics, 44(1), pp. 87-114, 2018.
% 
% Copyright (c) 2015-2019, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

v_ori = v;
TR = TriRep(f,v);
B = freeBoundary(TR);
bdy_index = B(:,1);

f_temp = f + length(v);

a = sort(bdy_index + length(v));
for i = length(a):-1:1
    f_temp(f_temp == a(i)) = a(i) - length(v);
end

f_new = [f;fliplr(f_temp)];
v_new = [v;v];

father = unique(f_new);
index = zeros(max(father),1);
index(father) = (1:size(father,1));
f_new = index(f_new);
v_ori_new = v_new(father,:);
v_new = v_ori_new;

inverse_father = zeros(max(father),1);
inverse_father(father) = (1:length(father))';
inverse_father(bdy_index+length(v)) = bdy_index;
L1 = cotangent_laplacian(v,f);
% combine the laplacian
[row,col,value] = find(L1);
newrow = [row; inverse_father(row+length(v))];
newcol = [col; inverse_father(col+length(v))];
newvalue = [value;value];
L = sparse(newrow, newcol, newvalue, length(v_new),  length(v_new));

vmap = sphericalmap_modified(v_new,f_new,v_new,L);

v_half = vmap(1:length(v),1:3);

%% using mobius transformation to place the face on the xy plane
% and then use stereographic projection and scaling

% map the boundary to be parallel to the xy plane
v_temp2 = v_half(bdy_index(2:end),:) - kron(ones(length(bdy_index)-1,1),v_half(bdy_index(1),:));
[~,index] = max(v_temp2(:,1).^2 + v_temp2(:,2).^2 +v_temp2(:,3).^2);

tri = v_half(bdy_index([1,index+1,round(index/2)+1]),:);

length1 = norm(tri(1,:) - tri(2,:),2);

x1 = 0; y1 = 0; x2 = length1; y2 = 0;

a = tri(2,1:3) - tri(1,1:3);
b = tri(3,1:3) - tri(1,1:3);

sin1 = (norm(cross(a,b),2))/(norm(a,2)*norm(b,2));
ori_h = norm(b,2)*sin1;
ratio = norm([x1-x2,y1-y2],2)/norm(a,2);
y3 = ori_h*ratio;
x3 = sqrt(norm(b,2)^2*ratio^2-y3^2);
x1 = x1 - length1/2; x2 = x2 -length1/2; x3 = x3 -length1/2;

plane1 = complex(v_half(:,1)./(1-v_half(:,3)), v_half(:,2)./(1-v_half(:,3)));
point1 = complex(v_half(bdy_index([1,index+1,round(index/2)+1]),1)./(1-v_half(bdy_index([1,index+1,round(index/2)+1]),3)), v_half(bdy_index([1,index+1,round(index/2)+1]),2)./(1-v_half(bdy_index([1,index+1,round(index/2)+1]),3)));

p1 = point1(1);
p2 = point1(2);
p3 = point1(3);

q1 = complex(x1,y1);
q2 = complex(x2,y2);
q3 = complex(x3,y3);

a = det([p1*q1,q1,1; p2*q2,q2,1; p3*q3,q3,1]);
b = det([p1*q1,p1,q1; p2*q2,p2,q2; p3*q3,p3,q3]);
c = det([p1,q1,1; p2,q2,1; p3,q3,1]);
d = det([p1*q1,p1,1; p2*q2,p2,1; p3*q3,p3,1]);

mobius = (a*plane1(:)+b)./(c*plane1(:)+d);
X = real(mobius(:));
Y = imag(mobius(:));

v_half_new = [2*X./(1+X.^2+Y.^2) , 2*Y./(1+X.^2+Y.^2) , (-1+X.^2+Y.^2)./(1+X.^2+Y.^2)];

% translate the boundary onto xy plane
v_half_new(:,3) = v_half_new(:,3) - mean(v_half_new(bdy_index,3));

% find a suitable stereographic projection
if mean(v_half_new(:,3))>=0
    v_disk = [v_half_new(:,1)./(1+v_half_new(:,3)), v_half_new(:,2)./(1+v_half_new(:,3)),0*X];
else
    v_disk = [v_half_new(:,1)./(1-v_half_new(:,3)), v_half_new(:,2)./(1-v_half_new(:,3)),0*X];
end

% suitable scaling and translating
v_disk = v_disk./(max(v_disk(:,1)) - min(v_disk(:,1)))*2;
v_disk(:,1) = v_disk(:,1) - (min(v_disk(:,1))+1);
v_disk(:,2) = v_disk(:,2) - (min(v_disk(:,2))+1);

%% Projection to unit circle
z_final = complex(v_disk(:,1),v_disk(:,2));
z_final(bdy_index) = z_final(bdy_index)./abs(z_final(bdy_index));
map = [real(z_final), imag(z_final), 1+0*z_final];

%% Correct distortion
mu = beltrami_coefficient(map, f, v_ori); 
map = linear_beltrami_solver(map,f,mu,bdy_index,map(bdy_index,1:2)); 
map(:,1) = -map(:,1);
end

