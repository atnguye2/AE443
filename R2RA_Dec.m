function [ alfa ,delta ] = R2RA_Dec( R );
delta = asind(R(3)/norm(R));                    % Declination
% Right ascension:
if ((R(2)/norm(R)) > 0)
    alfa = acosd((R(1)/norm(R))/cosd(delta));
else
    alfa = 360 - acosd((R(1)/norm(R))/cosd(delta));
end
end