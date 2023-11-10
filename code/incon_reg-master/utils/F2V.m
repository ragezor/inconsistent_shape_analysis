% -------------------------------------------------------
function S = F2V(v,f)
    ring = vertexAttachments(TriRep(f',v'));
    nv = length(v); nf = length(f);
    II = cellfun(@times,ring,num2cell(zeros(nv,1)),'UniformOutput',0);
    II = cell2mat(cellfun(@plus,II,num2cell(1:nv)','UniformOutput',0)')';
    JJ = cell2mat(ring')';
    avg = cellfun(@length,ring);
    S = sparse(II,JJ,ones(length(JJ),1),nv,nf);
    S = sparse(1:nv,1:nv,1./avg)*S;
end

