function b = in_or_out(bdy_xy,p)
% assume bdy_xy is oriented
[b,bb] = inpolygon(p(:,1),p(:,2),bdy_xy(:,1),bdy_xy(:,2));
b = b | bb;
end