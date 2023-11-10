% -------------------------------------------------------
function S = V2F(v,f)
    nv = length(v); nf = length(f);
    II = reshape(repmat(1:nf,3,1),3*nf,1);
    JJ = f(:);
    S = sparse(II,JJ,ones(length(JJ),1),nf,nv)./3;
end

