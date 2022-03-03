function fractures=periodicsplit(unsplitfracture,physdim,tol)
% Takes in a new fracture's vertices and DFN cell dimensions. Wraps the
% fractures around the cell to create a periodic geometry.

fractures=struct('points',[],'numpoints',0,'center',[],'normal',[],'size',0,'TriangleList',[]);

% Translate points (we want to work in only the positive quadrant)
points=unsplitfracture.points;
transpoints=points+repmat(physdim,size(points,1),1);

lx=physdim(1); ly=physdim(2); lz=physdim(3);


for i=1:27
    
    % set x range
    switch i
        case {1,4,7,10,13,16,19,22,25}
            xrange=[0 lx];
        case {2,5,8,11,14,17,20,23,26}
            xrange=[lx 2*lx];
        case {3,6,9,12,15,18,21,24,27}
            xrange=[2*lx 3*lx];
    end
    
    % set y range
    switch i
        case {1,2,3,10,11,12,19,20,21}
            yrange=[0 ly];
        case {4,5,6,13,14,15,22,23,24}
            yrange=[ly 2*ly];
        case {7,8,9,16,17,18,25,26,27}
            yrange=[2*ly 3*ly];
    end
    
    % set z range
    switch i
        case {1,2,3,4,5,6,7,8,9}
            zrange=[0 lz];
        case {10,11,12,13,14,15,16,17,18}
            zrange=[lz 2*lz];
        case {19,20,21,22,23,24,25,26,27}
            zrange=[2*lz 3*lz];
    end
    
    clippoints=cellclip(transpoints,xrange,yrange,zrange,tol);
    
    if ~isempty(clippoints)
        shift=min([xrange' yrange' zrange']);
        %         shift=physdim;
        fractures(end+1)=unsplitfracture;
        points=clippoints-repmat(shift,size(clippoints,1),1);
        points=arrangenodes(points,fractures(end).normal);
        fractures(end).points=points;
        numpoints=size(points,1);
        numtriangles=numpoints-2;
        trilist=[ones(numtriangles,1),(2:(numpoints-1))',(3:numpoints)'];
        fractures(end).TriangleList=trilist;
    end
    
end

fractures=fractures(2:end); % removed empty set

end



function clippoints=cellclip(points,xrange,yrange,zrange,tol)
xmin=xrange(1); xmax=xrange(2);
ymin=yrange(1); ymax=yrange(2);
zmin=zrange(1); zmax=zrange(2);

planepoint=[xmin 0 0]; planedirection=[1 0 0];
clippoints=polygonplaneclip(points,planepoint,planedirection,tol);

planepoint=[xmax 0 0]; planedirection=[-1 0 0];
clippoints=polygonplaneclip(clippoints,planepoint,planedirection,tol);

planepoint=[0 ymin 0]; planedirection=[0 1 0];
clippoints=polygonplaneclip(clippoints,planepoint,planedirection,tol);

planepoint=[0 ymax 0]; planedirection=[0 -1 0];
clippoints=polygonplaneclip(clippoints,planepoint,planedirection,tol);

planepoint=[0 0 zmin]; planedirection=[0 0 1];
clippoints=polygonplaneclip(clippoints,planepoint,planedirection,tol);

planepoint=[0 0 zmax]; planedirection=[0 0 -1];
clippoints=polygonplaneclip(clippoints,planepoint,planedirection,tol);

end