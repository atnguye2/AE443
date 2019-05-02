function alfa = Vincenty(lam1,phi1,lam2,phi2)
dellam = abs(lam1-lam2);
num = sqrt((cosd(phi2).*sind(dellam)).^2+(cosd(phi1).*sind(phi2)-sind(phi1).*cosd(phi2).*cosd(dellam)).^2);
den = sind(phi1).*sind(phi2)+cosd(phi1).*cosd(phi2).*cosd(dellam);
alfa = atan2d(num,den);
end