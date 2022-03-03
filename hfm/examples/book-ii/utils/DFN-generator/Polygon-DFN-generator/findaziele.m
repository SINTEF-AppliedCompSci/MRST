function [azimuth,elevation]=findaziele(direction,tol,varargin)

opt=struct('unit','radians');
opt=merge_options(opt,varargin{:});


% convert direction to column vector if it is a row
if size(direction,1)==1
    direction=direction';
end
assert(size(direction,1)==3);
assert(size(direction,2)==1);

% determine quadrant in x-y plane
xdir=[1;0;0];
ydir=[0;1;0];
zdir=[0;0;1];
xcoord=dot(direction,xdir);
ycoord=dot(direction,ydir);

azimuth=atan(ycoord/xcoord);

% correct for quadrant
quadrant=-1;
if xcoord<=tol && ycoord<=tol
    quadrant=0;
    azimuth=0; % prevent NaN output. In this case, elevation will be 0 or pi
elseif xcoord>=0 && ycoord>=0
    quadrant=1;
elseif xcoord<=0 && ycoord>=0
    quadrant=2;
    azimuth=pi+azimuth;
elseif xcoord<=0 && ycoord<=0
    quadrant=3;
    azimuth=pi+azimuth;
elseif xcoord>=0 && ycoord<=0
    quadrant=4;
    azimuth=2*pi+azimuth;
end

assert(quadrant>=0);


length=norm(direction);

zcoord=dot(direction,zdir);
elevation=acos(zcoord/length);

if strcmp(opt.unit,'degrees')
    azimuth=azimuth*180/pi;
    elevation=elevation*180/pi;
end



end