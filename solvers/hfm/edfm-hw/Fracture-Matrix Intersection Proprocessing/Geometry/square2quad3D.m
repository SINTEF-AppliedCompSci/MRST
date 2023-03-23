function quadpoints=square2quad3D(squarepoints,quadverts,tol,varargin)
% This function takes squarepoints, quadverts and transforms the square
% points into points in the quadrilateral defined by quadverts.

opt=struct('area',-1);
opt=merge_options(opt,varargin{:});

[maxdiag,mindiag,~,~]=findquaddiag3D(quadverts);

vertA=maxdiag(1); vertB=mindiag(1); vertC=maxdiag(2); vertD=mindiag(2);

%% scale square points to unit square
if opt.area<0
    area=convexpolygonarea(quadverts([vertA,vertB,vertC,vertD],:),tol);
    edgelength=realsqrt(area);
else
    edgelength=realsqrt(opt.area);
end

squarepoints=squarepoints./edgelength; %scaled to unit square

%% Get back quadrilateral points
W=quadverts(vertC,:)-quadverts(vertA,:);
U=quadverts(vertB,:)-quadverts(vertA,:);
V=quadverts(vertD,:)-quadverts(vertA,:);

UdotU=dot(U,U); VdotV=dot(V,V); UdotV=dot(U,V);
UdotW=dot(U,W); VdotW=dot(V,W);

A=[UdotU,UdotV;UdotV,VdotV]; b=[UdotW;VdotW];

sol=A\b;

a0=sol(1);
a1=sol(2);

numpoints=size(squarepoints,1);
quadpoints=zeros(numpoints,3);
for i=1:numpoints
    x0=squarepoints(i,1);
    x1=squarepoints(i,2);
    [y0,y1]=transform(x0,x1,a0,a1);
    quadpoints(i,:)=quadverts(vertA,:)+(y0*U)+(y1*V);   
end


end

function [y0,y1]=transform(x0,x1,a0,a1)
y0=(a0*x0)/((a0+a1-1)+(1-a1)*x0+(1-a0)*x1);
y1=(a1*x1)/((a0+a1-1)+(1-a1)*x0+(1-a0)*x1);
end