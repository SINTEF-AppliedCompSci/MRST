function [fracplanes,P32gen]=DFNfracplanesgenerator(fracplanes,fracinput,physdim,exclzonemult,tol,varargin)
% Generates fracplanes based on input. The position is randomized via
% Poisson process. The normal via Fisher distribution. The size via power
% law scaling. P32 fracture density is used to stop the fracture generation
% process. An exclusion zone method is used to prevent the fracture planes
% from being too close to each other. The domain is expanded about
% in volume to prevent boundary effects. Generated P32 is also available as
% output.

t1=clock;

opt=struct('displayexpdomain',false,'displayDFN',false,'circle',false);
opt=merge_options(opt,varargin{:});

if ~isempty(fracplanes)
    tempstore=fracplanes;
    clear fracplanes;
end

numpoints=fracinput.vertices;
P32=fracinput.P32;
normal=fracinput.normal.direction;
normalK=fracinput.normal.K;
minsize=fracinput.size.minsize;
maxsize=fracinput.size.maxsize;
fractaldim=fracinput.size.fractaldim;
perm=fracinput.perm;
poro=fracinput.poro;
aperture=fracinput.aperture;

%% Expand boundary
expdomain=physdim+maxsize*ones(1,3); % expanded domain
shift=maxsize*ones(1,3); % translation vector

%% Generate stopping criterion
cellvol=physdim(1)*physdim(2)*physdim(3);
targetarea=P32*cellvol; % this is the stopping criterion


%% Generate DFN in expanded domain
% % initialize first fracture
% newloc=expdomain.*rand(1,3); % random location
% newnormal=randdirection(normal,normalK,tol); % random normal
% newsize=randomsize_powerlaw(minsize,maxsize,fractaldim); % random size
% oldfracplanes=generateregularfracture(numpoints,newsize,newloc,newnormal); % generate fracplanes(1)
% oldfracplanes=createexclusionzone(oldfracplanes,exclzonemult); % generate exclusionzone for fracplanes(1)
% fracplanes=oldfracplanes; % initialize to get structure
% points=clip(oldfracplanes.points,physdim,shift,tol);
% if ~isempty(points),
%     count=1; 
%     fracplanes.points=points; 
% else
%     count=0; 
% end
% cumularea=convexpolygonarea(points,tol);

oldfracplanes=struct('points',[],'numpoints',0,'center',[],'normal',[],...
    'size',0,'TriangleList',[],'exclusionzone',[],'perm',0,'poro',0,'aperture',0,'area',0); % setting template
fracplanes=oldfracplanes;
cumularea=0;
discarded=0;
count=1;
while cumularea<targetarea
    % Generate new fracture
    newloc=expdomain.*rand(1,3);
    newnormal=randdirection(normal,normalK,tol);
    newsize=randomsize_powerlaw(minsize,maxsize,fractaldim);
    newfracture=generateregularfracture(numpoints,newsize,newloc,newnormal,'circle',opt.circle);
    
    % Check if new fracture intersects any old fracture exclusion zones
    for i=1:length(oldfracplanes) % first fracplane is useless
        if i==1, xsect=false; continue; end
        xsect=intersectexclusionzone(newfracture,oldfracplanes(i),tol);
        if xsect, break; end
    end
    if xsect, discarded=discarded+1; continue; end
    
    
    % if the above is successfully executed without going to the next
    % while loop, then generate exclusion zone. Note that exclusion zone is
    % for expanded domain.
    newfracture=createexclusionzone(newfracture,exclzonemult);
    newfracture.perm=perm;
    newfracture.poro=poro;
    newfracture.aperture=aperture;
    newfracture.area=0; % just to conform to template
%     newfracture
%     oldfracplanes
    oldfracplanes(end+1)=newfracture;   
    
    % clip with original domain. If points
    % is not empty, then this fracture has a portion inside the domain. In
    % this case, save the new fracture to fracplanes and save the clipped
    % points.
    points=clip(newfracture.points,physdim,shift,tol);
    if ~isempty(points)
        count=count+1;
        fracplanes(count)=newfracture;
        fracplanes(count).points=points;
        
        % Compute area within domain
        area=convexpolygonarea(points,tol);
        fracplanes(count).area=area;
        cumularea=cumularea+area;
    end  
end

fracplanes=fracplanes(2:end); % remove first one
fracplanes=rmfield(fracplanes,{'exclusionzone','TriangleList'});

% append previous fracplanes data
if exist('tempstore','var')
    fracplanes=[tempstore,fracplanes];
end

%% Visualization

if opt.displayexpdomain
    figure;
    for i=1:length(oldfracplanes)
        points=oldfracplanes(i).points;
        X=points(:,1);
        Y=points(:,2);
        Z=points(:,3);
        fill3(X,Y,Z,[1 1 0]);
        hold on;   
    end
    set(gca,'Zdir','reverse')
    axis equal tight;
    view(45,30);
    hold off;
end

if opt.displayDFN
    figure;
    for i=1:length(fracplanes)
        points=fracplanes(i).points;
        X=points(:,1);
        Y=points(:,2);
        Z=points(:,3);
        fill3(X,Y,Z,[1 1 0]);
        hold on;   
    end
    set(gca,'Zdir','reverse')
    axis equal tight;
    view(45,30);
    hold off;
end



%% Verbose
disp([num2str(discarded),' fractures were discarded due to intersection with exclusion zone of existing fractures.']);

P32gen=cumularea/cellvol;
disp(['Generated P32 in domain is ',num2str(P32gen)]);

t2=clock;
e=etime(t2,t1);
disp(['Generation of DFN took ',num2str(e),' seconds.']);
end

function points=clip(points,celldim,shift,tol)
% clipping plane and direction to keep
planepointleft=[0 0 0]+shift; planedirectionleft=[1 0 0];
planepointright=celldim+shift; planedirectionright=[-1 0 0];
planepointtop=[0 0 0]+shift; planedirectiontop=[0 0 1];
planepointbottom=celldim+shift; planedirectionbottom=[0 0 -1];
planepointfront=[0 0 0]+shift; planedirectionfront=[0 1 0];
planepointback=celldim+shift; planedirectionback=[0 -1 0];

% clip against each plane
points=polygonplaneclip(points,planepointleft,planedirectionleft,tol);
if size(points,1)==0,return; end;
points=polygonplaneclip(points,planepointright,planedirectionright,tol);
if size(points,1)==0,return; end;
points=polygonplaneclip(points,planepointtop,planedirectiontop,tol);
if size(points,1)==0,return; end;
points=polygonplaneclip(points,planepointbottom,planedirectionbottom,tol);
if size(points,1)==0,return; end;
points=polygonplaneclip(points,planepointfront,planedirectionfront,tol);
if size(points,1)==0,return; end;
points=polygonplaneclip(points,planepointback,planedirectionback,tol);
if size(points,1)==0,return; end;
% if we pass through all clipping steps above and did not go to next
% loop, then we clipped a polygon

points=points-repmat(shift,size(points,1),1); % translate back relative to origin (compensate for expanded domain)

end