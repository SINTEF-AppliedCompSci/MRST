function [fracplanes,P32gen]=periodicDFNfracplanesgenerator(fracplanes,fracinput,physdim,exclzonemult,tol,varargin)
% Generates fracplanes based on input. The position is randomized via
% Poisson process. The normal via Fisher distribution. The size via power
% law scaling. P32 fracture density is used to stop the fracture generation
% process. An exclusion zone method is used to prevent the fracture planes
% from being too close to each other. A periodic geometry is enforced by
% wrapping fracture around. Generated P32 is also available as output.

t1=clock;

opt=struct('displayDFN',false,'circle',false);
opt=merge_options(opt,varargin{:});

if isempty(fracplanes)
    newsetID=1;
    newtwinsetID=1;
else
    newsetID=max(vertcat(fracplanes.SetID))+1;
    newtwinsetID=max(vertcat(fracplanes.twinsetID))+1;
    tempstore=fracplanes;
    clear fracplanes;
end

numpoints=fracinput.vertices;
P32=fracinput.P32;
normal=fracinput.normal.direction;
normalK=fracinput.normal.K;
minsize=fracinput.size.minsize;
maxsize=fracinput.size.maxsize;
exponent=fracinput.size.exponent;
perm=fracinput.perm;
poro=fracinput.poro;
aperture=fracinput.aperture;



%% Generate stopping criterion
cellvol=physdim(1)*physdim(2)*physdim(3);
targetarea=P32*cellvol; % this is the stopping criterion


%% Generate DFN in expanded domain
% % initialize first fracture

fracplanes=struct('points',[],'numpoints',0,'center',[],'normal',[],...
    'size',0,'TriangleList',[],'exclusionzone',[],'perm',0,'poro',0,'aperture',0,'area',0,'SetID',newsetID,...
    'twinsetID',0); % setting template
cumularea=0;
discarded=0;
tries=0;

while cumularea<targetarea
    % Generate new fracture
    newloc=physdim.*rand(1,3);
    newnormal=randdirection(normal,normalK,tol);
    newsize=randomsize_powerlaw(2*minsize,2*maxsize,exponent)/2;
    newfracture=generateregularfracture(numpoints,newsize,newloc,newnormal,'circle',opt.circle);
    newfractures=periodicsplit(newfracture,physdim,tol);
    
    % Check if new fracture intersects any old fracture exclusion zones
    assert(tries<101,'Tries exceeded 100.');
    for i=1:length(fracplanes) % first fracplane is useless
        if i==1, xsect=false; continue; end
        for j=1:length(newfractures)
            xsect=intersectexclusionzone(newfractures(j),fracplanes(i),tol);
            if xsect, break; end
        end
        if xsect, break; end
    end
    if xsect, discarded=discarded+1; tries=tries+1; continue; end
    tries=0;
    
    % if the above is successfully executed without going to the next
    % while loop, then generate exclusion zone.
    for i=1:length(newfractures)
        newfractures(i).exclusionzone=[];
        newfractures(i)=createexclusionzone(newfractures(i),exclzonemult);
        newfractures(i).perm=perm;
        newfractures(i).poro=poro;
        newfractures(i).aperture=aperture;
        newfractures(i).area=convexpolygonarea(newfractures(i).points,tol); % just to conform to template
        newfractures(i).SetID=newsetID;
        cumularea=cumularea+newfractures(i).area;
    end
    
    if length(newfractures)>1
        [newfractures.twinsetID]=deal(newtwinsetID);
        newtwinsetID=newtwinsetID+1;
    else
        [newfractures.twinsetID]=deal(0);
    end

    fracplanes=[fracplanes,newfractures];
    
    disp(['Fracture #',num2str(length(fracplanes)-1),' successfully generated! ',num2str(round(100*cumularea/targetarea,1,'decimals')),'% completed.']);
 
end

fracplanes=fracplanes(2:end); % remove first one
fracplanes=rmfield(fracplanes,{'exclusionzone','TriangleList'});

% append previous fracplanes data
if exist('tempstore','var')
    fracplanes=[tempstore,fracplanes];
end

%% Visualization

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
%     xlim([0 physdim(1)]); ylim([0 physdim(2)]); zlim([0 physdim(3)]);
    plotGrid(cartGrid([1 1 1],physdim),'facealpha',0);
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

