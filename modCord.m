function [long,lat] = modCord(long,lat)
% long(or(lat>90,lat<-90)) = long(or(lat>90,lat<-90))+180;
% lat(lat>90) = 91-mod(lat(lat>90),90);
% lat(lat<-90) = -90-mod(lat(lat<-90),-90);

if lat > 90
    long(lat>90) = long(lat>90)+180;
    lat(lat>90) = 180-lat(lat>90);
end
if lat < -90
    long(lat<-90) = long(lat<-90)+180;
    lat(lat<-90) = -180-lat(lat<-90);
end
long = mod(long,360);
end