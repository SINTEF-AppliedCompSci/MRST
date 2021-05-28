function model = upscaleModelTPFA(model, partition, varargin)
%Upscale a fine model into a coarser version using a partition vector
%
% SYNOPSIS:
%   modelcoarse = upscaleModelTPFA(model, partition)
%
% DESCRIPTION:
%   Given a fine model and a partition vector (created using standard
%   coarsegrid tools such as partitionUI), this routine makes a coarser
%   model with upscaled properties. Using somewhat reasonable defaults,
%   most parts of the routine can be overriden by better values if the user
%   already knows for instance transmissibilities.
%
% REQUIRED PARAMETERS:
%   model     - Fine scale model. Subclass of `ReservoirModel`.
%
%   partition - Partition vector. Length equal to `model.G.cells.num`, with
%               positive indicator values. All cells with the same
%               indicator will be combined into a single coarse block.
%
% OPTIONAL PARAMETERS:
%   'validatePartition' - Ensure partition is connected on the grid, and
%                         numbered from 1...N without gaps.
% 
%   'transCoarse'       - Coarse transmissibilities. Will be calculated from
%                         upscaled permeability if not provided.
%
%   'transFromRock'     - Compute transmissibility from rock. Default is
%                         true. If disabled, the coarse transmissibility
%                         will be a sum instead.
%
%   'permCoarse'        - Coarse permeability. Will be calculated using
%                         harmonic averaging if not provided.
%
%   'neighborship'      - Coarse neighborship (matching transCoarse). Will
%                         be derived from fine grid if not provided.
%
%   'poroCoarse'        - Coarse porosities. Computed from fine model using
%                         a simple sum if not provided.
%
% RETURNS:
%   model  - Coarse model.
%
% SEE ALSO:
%   `upscaleSchedule`, `upscaleState`, `generateCoarseGrid`

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    opt = struct('validatePartition', true,...
                 'transFromRock',     true,...
                 'transCoarse',       [], ...
                 'permCoarse',        [], ...
                 'neighborship',      [], ...
                 'poroCoarse',        []);
     
    opt = merge_options(opt, varargin{:});
    
    % Handle grid
    require coarsegrid
    
    G = model.G;
    rock = model.rock;
    
    CG = getGrid(G, partition, opt);
    
    rock_c = getRock(rock, CG, opt);
    
    [Tc, Nc] = getTransmissibility(CG, rock_c, opt);
    if ~opt.transFromRock && isempty(opt.transCoarse)
        cts = rldecode((1:CG.faces.num)', diff(CG.faces.connPos));
        Tc = accumarray(cts, model.operators.T_all(CG.faces.fconn));
    end
    model.G = CG;
    model.rock = rock_c;
    model = model.setupOperators(CG, rock_c, 'neighbors', Nc, 'trans', Tc);
    model.FlowDiscretization = [];
end

function CG = getGrid(G, partition, opt)
    if isstruct(partition)
        CG = partition;
        return
    end
    if opt.validatePartition
        partition = processPartition(G, partition);
        partition = compressPartition(partition);
    end
    CG = generateCoarseGrid(G, partition);
    CG = coarsenGeometry(CG);
end

function rock_c = getRock(rock, CG, opt)
    % Handle rock
    poro_c = opt.poroCoarse;
    perm_c = opt.permCoarse;
    
    p = CG.partition;
    vol = CG.parent.cells.volumes;
    coarsevol = accumarray(p, vol);
    counts = accumarray(p, 1);
    if isfield(rock, 'poro') && isempty(poro_c)
        poro = rock.poro;
        % Include NTG directly into porosity for the time being...
        if isfield(rock, 'ntg')
            poro = poro.*rock.ntg;
        end
        % Sum up pore volumes - we want upscaled model to have same porosity.
        poro_c = accumarray(p, poro.*vol)./coarsevol;
    end
    
    if isfield(rock, 'perm') && isempty(perm_c)
        nK = size(rock.perm, 2);
        perm_c = zeros(CG.cells.num, nK);
        for i = 1:nK
            perm_c(:, i) = accumarray(p, 1./rock.perm(:, i));
        end
        perm_c = 1./bsxfun(@rdivide, perm_c, counts);
    end
    rock_c = makeRock(CG, perm_c, poro_c);
end

function [Tc, N_int] = getTransmissibility(CG, rock_c, opt)
    N = opt.neighborship;
    if isempty(N)
        N = CG.faces.neighbors;
    end
    
    % Handle transmissibility
    Tc = opt.transCoarse;
    nT = numel(Tc);
    
    % No. half faces
    nHf = size(CG.cells.faces, 1);
    % Number of faces
    nF  = size(N, 1);
    % No. interfaces
    intx = all(N ~= 0, 2);
    nIF = sum(intx);
    
    if nT == 0 
        % We have no upscaled transmissibility given we calculated them from the upscaled permeability    
        Tc = getFaceTransmissibility(CG, rock_c);
    elseif  nT == nF 
        % We have transmissibility in the format
        % one entry per face. This case is fine for setupOperators.
    elseif nT == nHf
        Tc  = 1 ./ accumarray(CG.cells.faces(:,1), 1./Tc, [nF, 1]);
    elseif nT == nIF
        % Use upscaled perm to compute trans for boundary and the provided
        % trans for internal faces.
        Tc_dummy = computeTrans(CG, rock_c);
        Tc_intf = Tc;
        
        Tc = zeros(nF, 1);
        Tc( intx) = Tc_intf(intx);
        Tc(~intx) = Tc_dummy(~intx);
    else
        msg = ['Number of input transmissibility entries (', num2str(nT), ...
            ') does not match either half face count (' num2str(nHf) '),', ...
            ' interface count (', num2str(nIF), ') or facecount (', num2str(nF), ').'];
        error(msg);
    end
    N_int = N(intx, :);
end
