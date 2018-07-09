function area=convexpolygonarea(nodes,tol, varargin)
% POLYGON AREA CALCULATOR - This function takes in the coordinates of the
% vertices of the 3D polygon in question (n-by-3 matrix, each row
% containing xyz coordinates) and calculates the area of the polygon.
%
% Second argument is a 'pn'/'pv' pair to specify if the vertices are
% arranged (clockwise/counterclockwise). The default is that the vertices
% are NOT arranged. 'pn' for this purpose is 'arranged' and 'pv' is true
% or false. We can also supply a normal direction using the 'planenormal'
% argument. This is recommended.

numvertices=size(nodes,1);

if numvertices<3 % if polygon is a line, point or has no vertices, area is 0
    area=0; return;
end

opt = struct('arranged',false,'normal',-1);
opt = merge_options(opt, varargin{:});
arranged=opt.arranged;

if opt.normal<-0.5
    normal=getnormal(nodes,tol);
else
    normal=opt.normal;
end

if ~arranged && numvertices>3 
    % arranges nodes
    nodes=arrangenodes(nodes,normal,'numvertices',numvertices);
end % if arranged, or is triangle, compute area directly
    
area=polygonareaarranged(nodes,numvertices,normal);

end

function area=polygonareaarranged(nodes,numvertices,normal)
% calculates area of an arranged polygon

% Van Gelder (1995)
h=floor((numvertices-1)/2);

if mod(numvertices,2)==0
    even=true;
else
    even=false;
end

if even
    Vk=nodes(numvertices-1,:);
else
    Vk=nodes(numvertices,:);
end

sum=0; % initialization
V0=nodes(numvertices,:);
V2h=nodes(2*h,:);
V2hminus1=nodes((2*h-1),:);
for i=1:(h-1)
    V2i=nodes(2*i,:);
    V2iplus1=nodes((2*i+1),:);
    V2iminus1=nodes((2*i-1),:);
    sum=sum+cross(V2i-V0,V2iplus1-V2iminus1);
end
sum=sum+cross(V2h-V0,Vk-V2hminus1);
area=abs(0.5*dot(normal,sum));

end








