function squarepoints=quad2square3D(quadpoints,tol,varargin)
% This function takes the vertices of a 3D quadrilateral and performs a
% perspective projection onto a square on the x-y plane. Area is preserved
% unless specified.
%
% Refer to https://www.geometrictools.com/Documentation/PerspectiveMappings.pdf

opt=struct('unitarea',false);
opt=merge_options(opt,varargin{:});

squarepoints=zeros(4,3);

[maxdiag,mindiag,~,~]=findquaddiag3D(quadpoints);

vertA=maxdiag(1); vertB=mindiag(1); vertC=maxdiag(2); vertD=mindiag(2);

if opt.unitarea
    edgelength=1;
else
    area=convexpolygonarea(quadpoints([vertA,vertB,vertC,vertD],:),tol);
    edgelength=realsqrt(area);
end

W=quadpoints(vertC,:)-quadpoints(vertA,:);
U=quadpoints(vertB,:)-quadpoints(vertA,:);
V=quadpoints(vertD,:)-quadpoints(vertA,:);

UdotU=dot(U,U); VdotV=dot(V,V); UdotV=dot(U,V);
UdotW=dot(U,W); VdotW=dot(V,W);

A=[UdotU,UdotV;UdotV,VdotV]; b=[UdotW;VdotW];

sol=A\b;

a0=sol(1); a1=sol(2);

%% vertA transformation
y0=0; y1=0;
[x0,x1]=transform(y0,y1,a0,a1,edgelength);
squarepoints(vertA,1:2)=[x0,x1];

%% vertB transformation
y0=1; y1=0;
[x0,x1]=transform(y0,y1,a0,a1,edgelength);
squarepoints(vertB,1:2)=[x0,x1];

%% vertC transformation
y0=a0; y1=a1;
[x0,x1]=transform(y0,y1,a0,a1,edgelength);
squarepoints(vertC,1:2)=[x0,x1];

%% vertD transformation
y0=0; y1=1;
[x0,x1]=transform(y0,y1,a0,a1,edgelength);
squarepoints(vertD,1:2)=[x0,x1];

end

function [x0,x1]=transform(y0,y1,a0,a1,edgelength)
x0=(a1*(a0+a1-1)*y0)/(a0*a1+a1*(a1-1)*y0+a0*(a0-1)*y1);
x1=(a0*(a0+a1-1)*y1)/(a0*a1+a1*(a1-1)*y0+a0*(a0-1)*y1);
x0=x0*edgelength;
x1=x1*edgelength;
end