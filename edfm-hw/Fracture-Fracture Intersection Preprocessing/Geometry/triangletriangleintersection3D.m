function [truth,line,length,df]=triangletriangleintersection3D(T1,T2,tol)
% This function determines if triangles intersect. It also returns, if 
% requested, the intersection line of two 3D triangular
% planes. Additionally, the length of the line is also returned. The last
% output is the subsegment centroid normal distances to the intersection
% line (See Moinfar et al 2013). The inputs are T1 and T2. These are 3x3
% matrices with each row being the x,y,z coordinates of the vertices of the
% triangles. Only works for non-parallel planes. If the lines do not
% intersect (truth==false),line=[], length=0, df=[]. Point
% intersection is considered as no intersection. line is a 2x3 matrix (each
% row is xyz coordinate of endpoint).

% Triangle 1
V1_1=T1(1,:);
V2_1=T1(2,:);
V3_1=T1(3,:);
normal_1=cross(V2_1-V1_1,V3_1-V1_1);
normal_1=normal_1/norm(normal_1);

% Triangle 2
V1_2=T2(1,:);
V2_2=T2(2,:);
V3_2=T2(3,:);
normal_2=cross(V2_2-V1_2,V3_2-V1_2);
normal_2=normal_2/norm(normal_2);

% check if T1 and T2 are parallel
paralleltest=norm(cross(normal_1,normal_2));
if paralleltest<tol
    % If triangles are parallel, they could overlap, not intersect, or
    % intersect at a point, or intersect on a line. For our case, we are
    % looking for full edge intersections. The specific case we deal with
    % are fracture segments of a curved fracture. So we will only look for
    % overlapping nodes
    Tall=[T1;T2];
    [~,IA]=uniquetol(Tall,tol,'ByRows',true,'DataScale',1);
    ixoverlaprows=setdiff(1:6,IA);
    xline=Tall(ixoverlaprows,:);
    
    if size(xline,1)~=2
        truth=false;line=[];length=0;df=[];
        return;
    end
else
    % Non parallel triangles
    [T1xP2,segmentT1xP2,typeT1xP2]=triangleplaneintersection3D(T1,normal_2,V1_2,tol);
    [T2xP1,segmentT2xP1,typeT2xP1]=triangleplaneintersection3D(T2,normal_1,V1_1,tol);
    
    % check if segments exist. If no segments, then no intersection line
    if T1xP2==false || T2xP1==false
        truth=false;line=[];length=0;df=[];
        return;
    end
    
    % check if any segment has zero length. Zero length means the intersection
    % between two triangles either does not exist, or is a point. Whichever
    % case, this is treated as no intersection.
    if typeT1xP2=='p' || typeT2xP1=='p'
        truth=false;line=[];length=0;df=[];
        return;
    end
    
    % since it is established that the segments exist, and that neither of them
    % are of zero length, we now test if the two segments intersect. In theory,
    % these lines should be co-linear, so we should not get an error from the
    % following function call.
    [T1xT2,xline]=truncatecolinearlines3D(segmentT1xP2,segmentT2xP1,tol);
    
    % check if the lines overlap
    if T1xT2==false
        truth=false;line=[];length=0;df=[];
        return;
    elseif (T1xT2==true) && (size(xline,1)==1) % overlap at one single point. Disregard.
        truth=false;line=[];length=0;df=[];
        return;
    end
end

truth=true;
line=xline;
length=norm(xline(1,:)-xline(2,:));
df=calctriangledf(T1,T2,xline,tol);



end


