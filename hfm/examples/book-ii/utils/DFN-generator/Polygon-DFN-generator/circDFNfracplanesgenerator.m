function fracplanes=circDFNfracplanesgenerator(numpoints,celldim,setdensity,setnormal,setnormalK, ...
    setradius,setradiusvar,shadowmult,tol)
% fracplane generator for EDFM. Function expands boundary by 30% in all
% directions and creates a larger space for DFN generation. Then circles
% are generated based on input parameters. Discrete points are generated to
% represent the circles, then clipped against the original domain
% dimensions. The output is a struct fracplanes which contains the field
% points. The field points will contain the vertices of a fracture plane in
% a sequential manner.

expdomain=celldim*1.6; % expanded domain
shift=celldim*0.3; % translation vector

% clipping plane and direction to keep
planepointleft=[0 0 0]+shift; planedirectionleft=[1 0 0];
planepointright=celldim+shift; planedirectionright=[-1 0 0];
planepointtop=[0 0 0]+shift; planedirectiontop=[0 0 1];
planepointbottom=celldim+shift; planedirectionbottom=[0 0 -1];
planepointfront=[0 0 0]+shift; planedirectionfront=[0 1 0];
planepointback=celldim+shift; planedirectionback=[0 -1 0];

circledfn=circleDFNgenerator(expdomain,setdensity,setnormal,setnormalK, ...
    setradius,setradiusvar,shadowmult,tol);

fracplanes=struct('points',[]);
count=0;
for i=1:length(circledfn)
    % generate polygon for each circle
    points=circlepoints(numpoints,circledfn(i).center,circledfn(i).normal,circledfn(i).radius,'randomize',true);
    
    % clip against each plane
    points=polygonplaneclip(points,planepointleft,planedirectionleft,tol);
    if size(points,1)==0,continue; end;
    points=polygonplaneclip(points,planepointright,planedirectionright,tol);
    if size(points,1)==0,continue; end;
    points=polygonplaneclip(points,planepointtop,planedirectiontop,tol);
    if size(points,1)==0,continue; end;
    points=polygonplaneclip(points,planepointbottom,planedirectionbottom,tol);
    if size(points,1)==0,continue; end;
    points=polygonplaneclip(points,planepointfront,planedirectionfront,tol);
    if size(points,1)==0,continue; end;
    points=polygonplaneclip(points,planepointback,planedirectionback,tol);
    if size(points,1)==0,continue; end;
    % if we pass through all clipping steps above and did not go to next
    % loop, then we clipped a polygon
    
    % save polygon in fracplanes
    count=count+1;
    points=points-repmat(shift,size(points,1),1); % translate back relative to origin (compensate for expanded domain)
    fracplanes(count).points=points;
end

end