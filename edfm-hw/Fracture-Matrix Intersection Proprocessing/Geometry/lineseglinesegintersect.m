function [truth,point] = lineseglinesegintersect(endpts1,endpts2,tol)
% Determines if two 3D line segments intersect.
% Output truth (true/false), point (one/two point or [])

dir1=endpts1(2,:)-endpts1(1,:);
dir2=endpts2(2,:)-endpts2(1,:);

% check if they are parallel. Parallel lines may either be collinear or do
% not intersect at all
unitdir1=dir1/norm(dir1);
unitdir2=dir2/norm(dir2);

if abs(abs(dot(unitdir1,unitdir2))-1)<tol %parallel case
    area1=heronsformula(endpts1(1,:),endpts1(2,:),endpts2(1,:),tol);
    area2=heronsformula(endpts1(1,:),endpts1(2,:),endpts2(2,:),tol);
    
    if area1<tol && area2<tol 
        %collinear case
        [overlap,point]=truncatecolinearlines(endpts1,endpts2,tol);
        if overlap %lines overlap
            truth=true;
        else
            truth=false;
        end        
    else
        %not collinear case
        truth=false; point=[];
    end
    return;
end

normal=cross(dir1,dir2);

normaldist=abs(dot(endpts2(1,:)-endpts1(1,:),normal));

% if the normal distance between two lines is non-zero, they do not
% intersect
if normaldist>tol
    truth=false; point=[]; return;
end

A=[dir1',-dir2'];
b=endpts2(1,:)'-endpts1(1,:)';
sol=A\b;

point=[endpts1(1,:)+dir1*sol(1);
    endpts2(1,:)+dir2*sol(2)];
point=uniquetol(point,tol,'ByRows',true,'DataScale',1);

% redundant check. Remove if desired.
if size(point,1)>1
    disp('Check that these two rows are the same. Ignore if same. If different, check function lineseglinesegintersect');
    disp(point);
    point=point(1,:);
end

if any(sol<-tol)||any(sol>(1+tol))
    truth=false;point=[]; return;
else
    truth=true;return;
end


end

