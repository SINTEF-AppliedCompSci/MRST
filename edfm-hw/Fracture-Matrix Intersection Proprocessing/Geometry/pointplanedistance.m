function dist=pointplanedistance(point,planenormal,planepoint)
% calculates the normal distance of a point to a plane

planenormal=planenormal/norm(planenormal); % make unit just in case
relposition=point-planepoint;

dist=abs(dot(relposition,planenormal));


end