function PointsInR = point2coverage(long,lat,r)
long=round(long);
lat=round(lat);
PointsInR = [];

%Right Half Sweep
coverageH = r;
SweepLongCoord = long;
SweepLatCoord = lat+r;

while coverageH > 0
    if InRadius
        newPoints = [[SweepLongCoord*ones(size(SweepLatCoord:-1:(SweepLatCoord-2*coverageH)))]',[SweepLatCoord:-1:(SweepLatCoord-2*coverageH)]'];
        PointsInR = [PointsInR;newPoints];
        SweepLongCoord = SweepLongCoord + 1;
    else
        SweepLatCoord = SweepLatCoord-1;
        coverageH = coverageH - 1;
    end
end

%Left Half Sweep
coverageH = r;
SweepLongCoord = long;
SweepLatCoord = lat+r;

while coverageH > 0
    if InRadius
        newPoints = [[SweepLongCoord*ones(size(SweepLatCoord:-1:(SweepLatCoord-2*coverageH)))]',[SweepLatCoord:-1:(SweepLatCoord-2*coverageH)]'];
        PointsInR = [PointsInR;newPoints];
        SweepLongCoord = SweepLongCoord - 1;
    else
        SweepLatCoord = SweepLatCoord-1;
        coverageH = coverageH - 1;
    end
end

    function b = InRadius
        b = pdist([SweepLongCoord,SweepLatCoord;long,lat]) <= r;
    end

end