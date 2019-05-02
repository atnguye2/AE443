function [alfa,delta] = GroundTrack(orbitAlt,RAAN,i,theta)
R_m = 3389.5;        % Mars' radius
mu = 4.282837e4;       % Mars' gravitational parameter [km^3/s^2]
J2 = 1960.45e-6;
we = 360/(24.6597*3600);      % Mars' rotation [deg/s]

% Satallites
rp    =  orbitAlt  + R_m;       % [km] Perigee Radius
ra    =  orbitAlt  + R_m;       % [km] Apogee Radius
% theta =  0;                 % [deg] True anomaly
% RAAN  =  60 ;          % [deg] Right ascension of the ascending node
% i     =  20 ;           % [deg] Inclination
omega =  0  ;         % [deg] Argument of perigee
a = (ra+rp)/2;               % Semimajor axis
e = (ra -rp)/(ra+rp) ;       % Eccentricity
h = (mu*rp*(1 + e))^0.5;     % Angular momentum
T = 2*pi*a^1.5/mu^0.5;       % Period
dRAAN = -(1.5*mu^0.5*J2*R_m^2/((1-e^2)*a^3.5))*cosd(i)*180/pi;
domega = dRAAN*(2.5*sind(i)^2 - 2)/cosd(i);

% Initial state
[R0, V0] = Orbital2State( h, i, RAAN, e,omega,theta);
[ alfa0 ,delta0 ] = R2RA_Dec( R0 );

ind = 1;
eps = 1E-9;
dt = 600;        % time step [sec]
ti = 0;

while(ti <= T)
    E = 2*atan(tand(theta/2)*((1-e)/(1+e))^0.5);
    M = E  - e*sin(E);
    t0 = M/(2*pi)*T;
    t = t0 + dt;
    M = 2*pi*t/T;
    E = keplerEq(M,e,eps);
    theta = 2*atan(tan(E/2)*((1+e)/(1-e))^0.5)*180/pi;
    RAAN  = RAAN  +  dRAAN*dt ;
    omega = omega + domega*dt;
    [R V] = Orbital2State( h, i, RAAN, e,omega,theta);
    % Considering Mars' rotation
    fi_mars = we*ti;
    Rot = [cosd(fi_mars), sind(fi_mars),0;...
        -sind(fi_mars),cosd(fi_mars),0;0,0,1];
    R = Rot*R;
    [ alfa(ind) ,delta(ind) ] = R2RA_Dec( R );
    ti = ti+dt;
    ind = ind + 1;
end
end