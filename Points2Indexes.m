function indxs = Points2Indexes(xgrid,ygrid,points)
for i = 1:size(points,1)
[points(i,1),points(i,2)] = modCord(points(i,1),points(i,2));
end
[~, indxs] = ismember( points, [xgrid(:) ygrid(:)], 'rows' );
if nnz(indxs) ~= numel(indxs)
    a=1+1
end
end