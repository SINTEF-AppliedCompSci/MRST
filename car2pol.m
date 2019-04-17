function [theta,r] = car2pol(x,y)
%convert (x,y) to (th, r)
%   theta belongs to [0 2*pi]

r=sqrt(x.^2+y.^2);
theta=(y>=0).*acos(x./r)+(y<0).*(2*pi-acos(x./r));
end

