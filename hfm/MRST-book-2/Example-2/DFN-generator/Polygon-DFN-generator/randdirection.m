function normal=randdirection(meannormal,K,tol)
% generates a random normal direction based on a mean normal direction and
% a fisher constant, K. The larger the K value, the smaller the variation
% in direction.

% find azimuth and elevation of mean normal
[azimuth,elevation]=findaziele(meannormal,tol);

angle1=fisherrandom(K); % equivalent to elevation
angle2=2*pi*rand; % equivalent to azimuth

% assuming mean normal is along [0 0 1]
normal=directionvector(angle2,angle1);

% rotate based on azimuth and elevation
normal=directionvector(azimuth,elevation,'refdir',normal);

end

