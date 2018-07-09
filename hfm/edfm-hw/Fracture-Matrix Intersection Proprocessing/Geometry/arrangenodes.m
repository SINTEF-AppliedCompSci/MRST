function nodes=arrangenodes(nodes,normal,varargin)

opt=struct('numvertices',0);
opt=merge_options(opt,varargin{:});

if opt.numvertices==0
    numvertices=size(nodes,1);
else
    numvertices=opt.numvertices;
end

dirx=nodes(2,:)-nodes(1,:);
dirx=dirx/norm(dirx); % make unit
diry=cross(normal,dirx);

center=sum(nodes)/numvertices;
relposition=nodes-repmat(center,numvertices,1);
xpos=dot(relposition,repmat(dirx,numvertices,1),2);
ypos=dot(relposition,repmat(diry,numvertices,1),2);

angle=zeros(numvertices,1);

for i=1:numvertices
    if xpos(i)>=0 && ypos(i)>=0  % quadrant 1
        angle(i)=atan(ypos(i)/xpos(i));
    elseif xpos(i)>=0 && ypos(i)<=0 % quadrant 4
        angle(i)=2*pi+atan(ypos(i)/xpos(i));
    else % quadrant 2&3
        angle(i)=pi+atan(ypos(i)/xpos(i));
    end % there shouldn't be a case of xpos(i)==0 since we use the center as reference. Division by 0 won't happen.
end

[~,index]=sort(angle);

nodes=nodes(index,:);
 
end