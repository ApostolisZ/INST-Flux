function dpyrdt = pyrode(t,y,MAL,PEP,CP,tt,v)

id = find(tt == t,1);
if isempty(id)
    ids = find(tt<t,1);
    id1 = ids(end);
    id2 = id1+1;
    MAL_MID = MAL(:,id1) + (t-tt(id1)).*(MAL(:,id2)-MAL(id1))./(tt(id2)-tt(id1));
    PEP_MID = PEP(:,id1) + (t-tt(id1)).*(PEP(:,id2)-PEP(id1))./(tt(id2)-tt(id1));
else
    MAL_MID = MAL(:,id);
    PEP_MID = PEP(:,id);
end

dpyrdt = (v(1)/CP).*PEP_MID + (v(2)/CP).*MAL_MID - ((v(3)+v(4)+v(5))/CP).*y;
end