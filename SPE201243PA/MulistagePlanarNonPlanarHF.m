function [fracplanes,frac_centroid_s]=MulistagePlanarNonPlanarHF(varargin)
%{
Create multi-stage hydraulic fractures geometry
%This function can create Non planar fracture of any orientation

!Assume all fractures in a stage have the same half-length and height
!Assume each hydraulic fracture plane is a rectangl
Author: Harun Rashid % 
Date modified 12/26/2019
%}

opt = struct('numStages'        ,  0,                      ...
             'fracSpacing'    , [],                      ...
             'numFracsPerStage'    , [],                      ...
             'fracHalfLength'   , [],                      ...
             'fracHeight'       , [],                      ...
             'theta',[],...
             'clusterSpacing'      , [],                      ...             
             'heelCoord'        , [10.0 20.0 15.0]);
opt = merge_options(opt, varargin{:});


%Handle constant input
if(isscalar(opt.fracSpacing))
    opt.fracSpacing=repmat(opt.fracSpacing,1,opt.numStages-1);
end
if(isscalar(opt.numFracsPerStage))
    opt.numFracsPerStage=repmat(opt.numFracsPerStage,1,opt.numStages);
end
if(isscalar(opt.fracHalfLength))
    opt.fracHalfLength=repmat(opt.fracHalfLength,1,opt.numStages);
end
if(isscalar(opt.fracHeight))
    opt.fracHeight=repmat(opt.fracHeight,1,opt.numStages);
end
if(isscalar(opt.clusterSpacing))
    opt.clusterSpacing=repmat(opt.clusterSpacing,1,opt.numStages);
end


if numel(opt.fracHalfLength)~= opt.numStages
    disp('[Error] Incompatible Input array size of fracHalfLength against numStages!');
end

%Create struct of fracture planes
fracplanes = struct;

frac_centroid=opt.heelCoord;
frac_centroid_s=[];
fracIdx=1;
for si=1:opt.numStages
    for fi=1:opt.numFracsPerStage(si)
        if(fi>1), frac_centroid(1) = frac_centroid(1) + opt.clusterSpacing(si); end
        fracplanes(fracIdx).points = ...
            createNONPlanarFrac(frac_centroid,...
                                 opt.theta,...
                                 opt.fracHalfLength(si),...
                                 opt.fracHeight(si));
        fracIdx=fracIdx+1;
        frac_centroid_s=[frac_centroid_s frac_centroid(1)];%HReDIT
    end
    if(si~=opt.numStages), frac_centroid(1)= frac_centroid(1) + opt.fracSpacing(si); end
end


end

function [polys] = createNONPlanarFrac(centroid,theta, HalfLength, Height) %HReDIT 
%Create a Inclined/vertical hydraulic fractures geometry based on theta
polys=[centroid(1)-(Height/2)*cos(theta), centroid(2)-HalfLength,centroid(3)+(Height/2)*sin(theta);
       centroid(1)-(Height/2)*cos(theta), centroid(2)+HalfLength,centroid(3)+(Height/2)*sin(theta);
       centroid(1)+(Height/2)*cos(theta), centroid(2)+HalfLength,centroid(3)-(Height/2)*sin(theta);
       centroid(1)+(Height/2)*cos(theta), centroid(2)-HalfLength,centroid(3)-(Height/2)*sin(theta);];
end
     