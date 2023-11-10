function map = Polygon2Square(v,f,corner)
% Disk-to-square map
%
% Input:
% v: nv x 3 vertex coordinates 
% f: nf x 3 triangulations
% corner : 4x1 vertex indices
% 4 - 3
% |   |
% 1 - 2
%
% Output:
% map: resultant mapping
%
% If you use this code in your own work, please cite the following paper:
% [1] T. W. Meng, G. P.-T. Choi and L. M. Lui, 
%     "TEMPO: Feature-Endowed Teichmüller Extremal Mappings of Point Clouds."
%     SIAM Journal on Imaging Sciences, 9(4), pp. 1922-1962, 2016.
%
% Copyright (c) 2015-2019, Gary Pui-Tung Choi
% https://scholar.harvard.edu/choi

TR = TriRep(f,v);
ring = vertexAttachments(TR);
B = freeBoundary(TR);
edge = edge_division(B,corner);
temp_edge = edge;
temp_corner = [corner;corner(1)];
for i = 1:4
    for j = 1:4
        if sum(ismember(temp_corner(i:i+1),edge{j})) == 4
            temp_edge{i} = edge{j};
        end
    end
end

vertical = [[edge{4};edge{2}],[zeros(length(edge{4}),1);ones(length(edge{2}),1)]];
horizontal = [[edge{1};edge{3}],[zeros(length(edge{1}),1);ones(length(edge{3}),1)]];


map = map_reconstruct(zeros(length(f),1),v,f,ring,vertical,horizontal);