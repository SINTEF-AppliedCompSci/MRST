function dist=pointplanedistance(point,planenormal,planepoint)
% calculates the normal distance of a point to a plane

planenormal=planenormal./realsqrt(sum(planenormal.^2,2)); % make unit just in case
relposition=point-planepoint;

dist=abs(dot(relposition',planenormal'))';


end