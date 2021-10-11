function [truth,point,type] = linesegmentplaneintersect(endpoints,planepoint,planenormal,tol)
% LINE SEGMENT VS PLANE INTERSECTION - Inputs are endpoints (2 rows
% each containing xyz coordinates of line endpoints), planepoint (row
% vector containing xyz coordinates of a point on plane) and planenormal
% (row vector containing xyz coefficient of normal vector). Returns a true
% or false statement. Also returns type of intersection: which can be
% 'none', 'interior', 'line', 'endpoint'. Lastly, returns the intersection
% point(s)(xyz coordinate row vector or [])

% Methodology: Determine the sides that the endpoints are on (relative to
% the plane point and plane normal). If one or both endpoints are on the
% plane, truth is true and type is 'endpoint'. If both signs are the same,
% truth is false and type is 'none'. If both signs are different, then truth
% is true and type is 'interior'.

% sign calculation
relposition=endpoints-repmat(planepoint,2,1);
sign=dot(relposition,repmat(planenormal,2,1),2);
sign(abs(sign)>0)=sign(abs(sign)>0)./abs(sign(abs(sign)>0));

% determine truth and type
if any(abs(sign)<tol)
    truth=true; point=endpoints(abs(sign)<tol,:); 
    if size(point,1)==1
        type='endpoint';
    else
        type='line';
    end
    return;
elseif abs(sum(sign))<tol
    truth=true; type='interior';
else
    truth=false; type='none'; point=[]; return;
end

% if successfully exit if-elseif-else block, find the intersection point
r=dot(planenormal,planepoint-endpoints(1,:))/dot(planenormal,endpoints(2,:)-endpoints(1,:));
point=endpoints(1,:)+r*(endpoints(2,:)-endpoints(1,:));

end

