function df=calcPEBIdf(P1,P2,xline,tol)
% This function calculates df for P1 and P2 as described in Moinfar et al
% (2013). The inputs P1 and P2 are Nx3 matrices in which each row
% represents the xyz coordinates of the vertices of 3D PEBI planes. xline is
% a 2x3 matrix in which each row represents the xyz coordinates of the
% endpoints of the 3D intersection line between the two PEBI planes. The
% intersection line must have a non-zero length.
%
% The output df will be a 2x1 vector with the first element containing the
% value of df for P1, and second element containing value of df for P2.

% Intersection line unit normal
normalx=xline(1,:)-xline(2,:);
normalx=normalx/norm(normalx);

% PEBI 1 unit normal
normal1=cross(P1(1,:)-P1(2,:),P1(1,:)-P1(3,:));
normal1=normal1/norm(normal1);

% PEBI 2 unit normal
normal2=cross(P2(1,:)-P2(2,:),P2(1,:)-P2(3,:));
normal2=normal2/norm(normal2);

normalP=[normal1;normal2];
P=cell(2,1);
P{1}=P1;P{2}=P2;

df=zeros(2,1);
for i=1:2
    % extract relevant unit normal and vertices
    normal=normalP(i,:);
    V=P{i};
    numverts=size(V,1);
    
    % find unit vector perpendicular to both plane normal and line normal
    % and call it refdir (this is an in-plane vector)
    refdir=cross(normal,normalx);
    
    % determine the side the vertices are on. sign is a Nx1 vector with
    % each row containing a real number. If this number is positive, the
    % corresponding vertex in V is on the positive side of the unit normal,
    % vice versa. We do not consider vertices on the intersection line
    % because this is not necessary; the two endpoints of the intersection
    % line are enough to describe that edge of the subsegment.
    sign=dot(V-repmat(xline(1,:),numverts,1),repmat(refdir,numverts,1),2);
    
    positiveside=V(sign>tol,:); positiveside=[positiveside;xline];
    negativeside=V(sign<-tol,:); negativeside=[negativeside;xline];
    
    positivecentroid=convexpolygoncentroid(positiveside,tol);
    negativecentroid=convexpolygoncentroid(negativeside,tol);
    
    positivenormaldistance=abs(dot(positivecentroid-xline(1,:),refdir));
    negativenormaldistance=abs(dot(negativecentroid-xline(1,:),refdir));
    
    df(i,1)=(positivenormaldistance+negativenormaldistance)/2;
end






end