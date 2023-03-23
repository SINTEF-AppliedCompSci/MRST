function [fracplanes,frac_centroid_s]=createMultiStageHFsCGEOR(varargin)
%{
Create multi-stage hydraulic fractures geometry

!Assume all fractures in a stage have the same half-length and height
!Assume each hydraulic fracture plane is a rectangle
!Assume all hydraulic fractures are vertical and orthogonal to the X axis

Arguments
---------
fl        --  NFx4 array of fracture segment length 
xy_wells  --  NFx2 array of center well location

Author:Bin Wang
Date: Nov.21.2018
%}
opt = struct('numStages'        ,  0,                      ...
             'fracSpacing'    , [],                      ...
             'numFracsPerStage'    , [],                      ...
             'fracHalfLength'   , [],                      ...
             'fracHeight'       , [],                      ...
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

%TO VERIFY: Why are you comparing fracHalfLength  with number of frac
%stages?
%Bin: This is used for the case of variable paramters input mode 
%Check input paramters
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
            createPlanarFrac(frac_centroid,...
                             opt.fracHalfLength(si),...
                             opt.fracHeight(si));
        fracIdx=fracIdx+1;
        frac_centroid_s=[frac_centroid_s frac_centroid(1)];
    end
    if(si~=opt.numStages), frac_centroid(1)= frac_centroid(1) + opt.fracSpacing(si); end
end


end

function [polys] = createPlanarFrac(centroid, HalfLength, Height)
%Create a rectangle hydraulic fractures geometry
polys=[centroid(1), centroid(2)-HalfLength,centroid(3)-Height/2;
       centroid(1), centroid(2)+HalfLength,centroid(3)-Height/2;
       centroid(1), centroid(2)+HalfLength,centroid(3)+Height/2;
       centroid(1), centroid(2)-HalfLength,centroid(3)+Height/2;];
end

     