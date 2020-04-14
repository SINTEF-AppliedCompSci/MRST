function [subforces, mappings] = getSubForces(forces, mappings)
    [subW, subBC, subSRC, keepW, keepBC, keepSRC] = deal([]); %#ok
    % Wells
    if ~isempty(forces.W)
        mappings.wells = struct();
        [subW, keepW]  = getSubWells(forces, mappings.cells);
    end
    mappings.wells.keep = keepW;
    % Boundary conditions
    if ~isempty(forces.bc)
        [subBC, keepBC] = getSubBC(forces, mappings.faces); %#ok
    end
    % Source terms
    if ~isempty(forces.src)
        [subSRC, keepSRC] = getSubSources(forces, mappings); %#ok
    end
    subforces = forces;
    subforces.W   = subW;
    subforces.bc  = subBC;
    subforces.src = subSRC;
end

%-------------------------------------------------------------------------%
function [subW, keep] = getSubWells(forces, map)
    subW = forces.W;
    keep = false(numel(subW), 1);    
    for i = 1:numel(subW)
        keepw = map.keep(subW(i).cells);
        if any(keepw)
            assert(all(keepw));
            subW(i).cells = map.renum(subW(i).cells);
            keep(i) = true;
        end
    end
    subW = subW(keep);
end

%-------------------------------------------------------------------------%
function [subBC, keep] = getSubBC(forces, map)
    subBC = forces.bc;
    keep = map.keep(subBC.face);
    if ~any(keep)
        subBC = [];
        return
    end
    subBC.face  = map.renum(subBC.face(keep));
    subBC.type  = subBC.type(keep);
    subBC.value = subBC.value(keep);
    subBC.sat   = subBC.sat(keep,:);
end

%-------------------------------------------------------------------------%
function getSubSources(varargin)
    error('Source terms not implemented yet')
end