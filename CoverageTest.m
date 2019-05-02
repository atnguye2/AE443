function h = CoverageTest
h = figure;
hold on
axis([0 360 -90 90]);
ylabel('Latitude [deg]');
xlabel('Longitude [deg]');
title('Projection of Circular Coverage Radius onto Long/Lat');
[c1x,c1y] = circle(180,0,30);
[c2x,c2y] = circle(0,0,30);
[c2x,c2y] = modCord(c2x,c2y);
[c3x,c3y] = circle(200,90,30);
[c3x,c3y] = modCord(c3x,c3y);
[c4x,c4y] = circle(180,-90,30);
[c4x,c4y] = modCord(c4x,c4y);
[c5x,c5y] = circle(270,-105,30);
[c5x,c5y] = modCord(c5x,c5y);
plot(c1x,c1y,'b.')
plot(c2x,c2y,'k.')
plot(c3x,c3y,'r.')
plot(c4x,c4y,'c.')
plot(c5x,c5y,'m.')
end