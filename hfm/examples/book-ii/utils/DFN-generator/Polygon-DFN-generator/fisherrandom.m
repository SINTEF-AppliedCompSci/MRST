function angle=fisherrandom(K,varargin)
% generates a random number between 0 and pi/2 based on fisher
% distribution. pn/pv pair can be input to change output to an angle in
% degrees instead of the default in radians. 'unit'/'degrees'/'radians'

opt=struct('unit','radians');
opt=merge_options(opt,varargin{:});

angle=acos(1+log(1-rand)/K);

if strcmp(opt.unit,'degrees')
    angle=angle*180/pi;
end


end