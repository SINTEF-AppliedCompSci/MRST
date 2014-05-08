function mpos = mapAxes(pos, MA)
% Transform lateral position from map to local coordinate system 
%
% SYNOPSIS:
%   p = mapAxis(pos, ma)
%
% PARAMETERS:
%   pos - nx2 vector of map coordinates that are to be transformed
%
%   ma  - 6x1 vector of points: x1 y1 x2 y2 x3 y3
%         Here: (x1,y1) describes a point on the y-axis, (x2,y2) is the
%         grid origin, and (x3,y3) is a point on the x-axis
%
% RETURNS:
%   p   - nx2 vector of coordinates transformed to local coordinate system
%   
% EXAMPLE:
%   p = mapAxes([x y], [100 100 0 100 100 0])

%{
#COPYRIGHT#
%}

% $Date: 2012-12-20 12:29:31 +0100 (Thu, 20 Dec 2012) $
% $Revision: 10515 $

e1 = MA(5:6) - MA(3:4); e1 = e1/norm(e1);
e2 = MA(1:2) - MA(3:4); e2 = e2/norm(e2);
msign=dot([0 0 1],cross([e1,0],[e2,0]));
if(msign<0)
    warning('mapAxis','mapAxes may change signature of coordinate system.')
end
% tranform coordinate system may change the sign of the system
mpos = bsxfun(@times,pos(:,1),e1)+bsxfun(@times,pos(:,2),e2);
mpos = bsxfun(@plus,mpos,MA(3:4));
