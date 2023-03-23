function df=calctriangledf(T1,T2,xline,tol)
% This function calculates df for T1 and T2 as described in Moinfar et al
% (2013). The inputs T1 and T2 are 3x3 matrices in which each row
% represents the xyz coordinates of the vertices of 3D triangles. xline is
% a 2x3 matrix in which each row represents the xyz coordinates of the
% endpoints of the 3D intersection line between the two triangles. The
% intersection line must have a non-zero length.
%
% The output df will be a 2x1 vector with the first element containing the
% value of df for T1, and second element containing value of df for T2.
%
% Later improvement: code can be written such that T1 and T2's normals may
% be precomputed and used as inputs. This makes the overall preprocessing
% more efficient.

% Intersection line unit normal
normalx=xline(1,:)-xline(2,:);
normalx=normalx/norm(normalx);

% Triangle 1 unit normal
normal1=cross(T1(1,:)-T1(2,:),T1(1,:)-T1(3,:));
normal1=normal1/norm(normal1);

% Triangle 2 unit normal
normal2=cross(T2(1,:)-T2(2,:),T2(1,:)-T2(3,:));
normal2=normal2/norm(normal2);

normalT=[normal1;normal2];
T=cell(2,1);
T{1}=T1;T{2}=T2;

df=zeros(2,1);
for i=1:2
    % extract relevant unit normal and vertices
    normal=normalT(i,:);
    V1=T{i}(1,:);
    V2=T{i}(2,:);
    V3=T{i}(3,:);
    V=[V1;V2;V3];
    
    % find unit vector perpendicular to both plane normal and line normal
    % and call it refdir (this is an in-plane vector)
    refdir=cross(normal,normalx);
    
    % determine the side the vertices are on. sign is a 3x1 vector with
    % each row containing a real number. If this number is positive, the
    % corresponding vertex in V is on the positive side of the unit normal,
    % vice versa. We do not consider vertices on the intersection line
    % because this is not necessary; the two endpoints of the intersection
    % line are enough to describe that edge of the subsegment.
    sign=dot(V-repmat(xline(1,:),3,1),repmat(refdir,3,1),2);
    
    positiveside=V(sign>tol,:);
    negativeside=V(sign<-tol,:);
    
    positivecentroid=findtrisubsegmentcentroid(positiveside,xline);
    negativecentroid=findtrisubsegmentcentroid(negativeside,xline);
    
    positivenormaldistance=abs(dot(positivecentroid-xline(1,:),refdir));
    negativenormaldistance=abs(dot(negativecentroid-xline(1,:),refdir));
    
    df(i,1)=(positivenormaldistance+negativenormaldistance)/2;
end


end



