function direction=directionvector(azimuth,elevation,varargin)
% generates a direction vector from azimuth and elevation. pn/pv pair
% available to specify 'units' in 'degrees' or to change the reference
% direction 'refdir' to any other direction.
opt=struct('units','radians','refdir',[0 0 1]);
opt=merge_options(opt,varargin{:});

if strcmp(opt.units,'degrees')
    azimuth=azimuth*pi/180;
    elevation=elevation*pi/180;
end

direction=(opt.refdir)'; % make column vector for computation purposes
direction=rotaroundy(elevation)*direction;
direction=rotaroundz(azimuth)*direction;
direction=direction'; % make row vector per convention
end

function rotmat=rotaroundy(angle)
rotmat=[cos(angle) 0 sin(angle);
        0 1 0;
        (-sin(angle)) 0 cos(angle)];
end

function rotmat=rotaroundz(angle)
rotmat=[cos(angle) (-sin(angle)) 0;
        sin(angle) cos(angle) 0;
        0 0 1];
end