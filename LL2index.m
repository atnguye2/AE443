function ind = LL2index(long,lat,xgrid,ygrid)
[long,lat] = modCord(long,lat);
centerX = find(xgrid==round(long));
centerY = find(ygrid==round(lat));
ind = intersect(centerX,centerY);
end