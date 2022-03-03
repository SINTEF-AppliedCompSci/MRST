function truth=intersectexclusionzone(newfracture,oldfracture,tol)
% This function takes in the new and old fractures and check if the new
% fracture intrudes the exclusion zone of the old fracture. The output is a
% simple true/false statement.
%
% This function uses repeated triangle-triangle intersection checks to see
% if the new fracture intersects the old fracture's exclusion zone.

%% Quick filter by sphere
% Both fractures will be made into spheres. If the spheres do not collide,
% the fractures do not intersect.

newcenter=mean(newfracture.points,1);
pointrelpos=bsxfun(@minus,newfracture.points,newcenter);
newrad=realsqrt(max(sum(pointrelpos.^2,2)));

oldcenter=mean(oldfracture.points,1);
pointrelpos=bsxfun(@minus,oldfracture.exclusionzone.surface1.points,oldcenter);
oldrad=realsqrt(max(sum(pointrelpos.^2,2)));

centerdist=norm(newcenter-oldcenter);

% if the distance between centers is larger than the sum of radii, then the
% two spheres don't intersect, hence the two fractures don't intersect.
if centerdist > (newrad + oldrad)
    truth=false;
    return;
end



%% Check if any vertex of new fracture is within exclusion zone
numexzonesurf=length(fieldnames(oldfracture.exclusionzone)); % number of exclusion zone surfaces
newverts=newfracture.points;
numnewverts=size(newverts,1);
flag=zeros(numnewverts,1);
for i=1:numexzonesurf
    fieldname=['surface',num2str(i)];
    refpoint=oldfracture.exclusionzone.(fieldname).points(1,:);
    insidedir=oldfracture.exclusionzone.(fieldname).insidedirection;
    
    normdist=dot(newverts-repmat(refpoint,numnewverts,1),repmat(insidedir,numnewverts,1),2);
    normdist(normdist>0)=1;
    normdist(normdist<0)=-1;
    normdist(normdist==0)=0;
    
    flag=flag+normdist;
end

if any(flag==numexzonesurf)
    truth=true; return;
end



%% Exclusion zone surface intersection detection
% collect all triangles for new fracture regular polygon
trilist_new=newfracture.TriangleList;
points_new=newfracture.points;
numtri_new=size(trilist_new,1);
trinodes_new=cell(numtri_new,1);
for i=1:numtri_new
    trinodes_new{i}=points_new(trilist_new(i,:),:);
end

% collect all triangles for all surfaces of exclusion zone
numpoints_old=size(oldfracture.points,1);
numtri_old=4*(numpoints_old-1); % this can be easily worked out if we understand how the triangulation is done
trinodes_old=cell(numtri_old,1);
count=0;
for i=1:numexzonesurf
    fieldname=['surface',num2str(i)];
    trilist_old_i=oldfracture.exclusionzone.(fieldname).TriangleList;
    surfpoints_old_i=oldfracture.exclusionzone.(fieldname).points;
    numtri_surf_i=size(trilist_old_i,1);
    for j=1:numtri_surf_i
        count=count+1;
        trinodes_old{count}=surfpoints_old_i(trilist_old_i(j,:),:);
    end
end

% outer loop: exclusion zone
% inner loop: new fracture
for i=1:numtri_old
    T_old=trinodes_old{i};
    
    for j=1:numtri_new
        T_new=trinodes_new{j};
        
        %check for intersection
        [xsect,~,~]=triangletriangleintersection(T_new,T_old,tol);
        
        if xsect;
            truth=true; return;
        end
    end
end

truth=false;

end