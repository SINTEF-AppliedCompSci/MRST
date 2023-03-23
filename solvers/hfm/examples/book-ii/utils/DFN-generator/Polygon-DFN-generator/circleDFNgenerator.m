function circledfn=circleDFNgenerator(celldim,density,normal,normalK,radius,radiusvar,shadowmult,tol)
% DFN generator for one set of circular fractures. 
% This function takes in celldim which is the rectangular domain
% dimensions [xlength, ylength, zlength]. Density is the
% P31 density of the fracture set. normal is the direction of the mean
% normal direction. normalvar the variation of the normal direction.
% radius is the mean radius. radius var is the variation of the radius.
%
% circledfn (output) is a structure containing data on the radius, normal
% direction and center location of each fracture within celldim.
%
% Method used here is from Discontinuity Analysis for Rock Engineering by
% Priest (1993). Section 6.8.

% AS A START, THIS FUNCTION DOES NOT INCLUDE VARIATION IN NORMAL DIRECTION
% AND CIRCLE RADIUS.

cellvol=celldim(1)*celldim(2)*celldim(3);

targetnum=round(density*cellvol); % round to the nearest integer
fraccount=0;
circledfn=struct('center',[0 0 0],'normal',[0 0 0],'radius',0);
newcircle=struct('center',[0 0 0],'normal',[0 0 0],'radius',0);
intersectnum=0;
while fraccount<targetnum % stops when fraccount=targetnum
    [newcircle.center,newcircle.normal,newcircle.radius]=generateonecircle(celldim,normal,normalK,radius,radiusvar);
        
    % check for circle-shadowzone intersect with all previous circles
    for i=1:fraccount
        intersect=circleshadowintersect(newcircle,circledfn(i),shadowmult,tol);
        if (intersect), intersectnum=intersectnum+1; break; end
    end
    if (fraccount>0) && intersect, continue; end;
    
    % quick hack to add angle variation
    newcircle.normal=randdirection(normal,normalK,tol);
    
    
    fraccount=fraccount+1;
    circledfn(fraccount)=newcircle;
end

printout=['DFN completed with ',num2str(intersectnum),' fracture(s) discarded due to intersection with shadow zone.'];
disp(printout);

end

function [center,normal,radius]=generateonecircle(celldim,normal,normalvar,radius,radiusvar)
% generates a single circle
center=celldim.*rand(1,3); % random center location

end