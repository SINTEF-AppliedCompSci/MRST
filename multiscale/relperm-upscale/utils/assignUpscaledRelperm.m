function [model, kr] = assignUpscaledRelperm(model, kr, varargin)
%Assign upscaled relperm to AD model
%
% SYNOPSIS:
%   u[model, kr] = assignUpscaledRelperm(model, kr, varargin)
%
% REQUIRED PARAMETERS:
%   model    - Upscaled coarse model
%   kr       - Relperm struct (e.g. computeRelpermFromStates)      
%
%
% RETURNS:
%   model - Model with modified fluid object.
%   kr    - Regularized relperm object (equal sample points and unphysical
%           values removed.)
% SEE ALSO:
%   `regularizeSaturationFunction`, `computeRelpermFromStates`,
%   `mergeHalfFaceRelPerm`

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('resamplePoints', 25,... % Either explicit point list or number of points to resample
                 'minPoints', 3, ... % Minimum number of points for a relperm to be included
                 'setWells', true, ... % Set well relperm
                 'minCoverage', 0.0); % Points must cover this fraction of the [0, 1] saturation interval
    [opt, regopt] = merge_options(opt, varargin{:});
    
    numPh = numel(kr);
    for i = 1:numPh
        [model, kr{i}] = assignPhaseRelPerm(model, kr{i}, i, opt, regopt);
    end
end

function [model, relPermStruct]= assignPhaseRelPerm(model, relPermStruct, index, opt, regopt)
    nc = model.G.cells.num;
    kr = relPermStruct.reservoir;
    
    n_kr = numel(kr);
    nf = size(model.operators.N, 1);
    sz_kr = size(kr);
    
    % Set up uniform grid for re-sampling interpolants
    if isscalar(opt.resamplePoints)
        if isfinite(opt.resamplePoints)
            S = linspace(0, 1, opt.resamplePoints)';
        else
            % VERY EXPENSIVE
            all_S = cellfun(@(x) x.S, kr, 'UniformOutput', false);
            S = vertcat(all_S{:});
            S = [0; S; 1];
            S = unique(S);
        end
    else
        S = reshape(opt.resamplePoints, [], 1);
    end
    ns = numel(S);
    
    % Set up fluid relperm
    [phases, phaseNames] = model.getPhaseNames();
    
    kr_str = upper(phases(index));
    if strcmpi(kr_str, 'O')
        isOil = true;
        % Special case
        if isfield(model.fluid, 'krO')
            kr_str = 'O';
        elseif model.water
            kr_str = 'OW';
        else
            kr_str = 'OG';
        end
    else
        isOil = false;
    end
    kr_name = ['kr', kr_str];
    fn = model.fluid.(kr_name);
    
    
    count = 0;
    for i = 1:n_kr
        [kr{i}, bad] = processMissing(kr{i}, fn, opt, regopt);
        kr{i} = makeUniform(S, kr{i}, opt);
        count = count + double(bad);
    end
    krw = relPermStruct.wells;
    for i = 1:numel(relPermStruct.wells)
        for k = 1:numel(relPermStruct.wells{i}.perf)
            krw{i}.perf{k} = processMissing(krw{i}.perf{k}, fn, opt, regopt);
            krw{i}.perf{k} = makeUniform(S, krw{i}.perf{k}, opt);
        end
    end
    fprintf('Processing of phase %s OK. Replaced %d of %d functions due to insufficient data.\n', ...
            phaseNames{index}, count, numel(kr));
    
    tmp = cellfun(@(x) x.kr', kr, 'UniformOutput', false);
    f = [];
    for i = 1:size(tmp, 2)
        f = [f; vertcat(tmp{:, i})]; %#ok
    end
    if sz_kr(1) == nc
        % Cell-wise values are given
        reg = (1:nc)';
        relPerm = getUniformInterpRegLinear(S', f, reg);
    elseif sz_kr(1) == nf 
        offset = sz_kr(2)*nf;
        cell_region = repmat(offset+1, 1, nc);
        % Cell relperm is base relperm
        f = [f; fn(S)'];
        
        if sz_kr(2) > 1
            % Half-face values
            % First region for each half face, then for all cells
            reg = [(1:offset), cell_region]';
        else
            % Face values repeated twice for upwinding, then all cells
            reg = [(1:nf), (1:nf), cell_region]';
        end
        if opt.setWells
            wcells = [];
            wfcn = [];
            for wno = 1:numel(krw)
                fw = cellfun(@(x) x.kr', krw{wno}.perf, 'UniformOutput', false);
                fw = vertcat(fw{:});
                wc = krw{wno}.cells;
                
                wfcn = [wfcn; fw];
                wcells = [wcells; wc];
            end

            reg(offset+wcells) = (offset+1)+(1:numel(wcells));
            f = [f; wfcn];
        end
        relPermFace = getUniformInterpRegLinear(S', f, reg);
        
        N = model.operators.N;
        n1 = N(:, 1);
        n2 = N(:, 2);

        map = @(S) [S(n1); S(n2); S];
        relPerm = @(S, varargin) relPermFace(map(S));
    else
        error('Unable to map');
    end
    
    
    if isOil
        names = {'krOW', 'krO', 'krOG'};
        for i = 1:numel(names)
            if isfield(model.fluid, names{i})
                model.fluid.(names{i}) = relPerm;
            end
        end
    else
        model.fluid.(kr_name) = relPerm;
    end
    relPermStruct.reservoir = kr;
    relPermStruct.wells = krw;
end

function [kr, bad] = processMissing(kr, replacefn, opt, regopt)
    ok = checkCoverage(kr, opt);
    
    [kr.S, kr.kr] = regularizeSaturationFunction(kr.S, kr.kr, regopt{:});
    % Re-check after processing
    ok = ok && checkCoverage(kr, opt);
    bad = ~ok;
    if bad
        kr.S = (0:0.01:1)';
        kr.kr = replacefn(kr.S);
    end
end

function ok = checkCoverage(kr, opt)
    ok = numel(kr.S) >= opt.minPoints &&... % Too few points given
          abs(kr.S(end) - kr.S(1)) >= opt.minCoverage; % Too small interval covered

end

function kr = makeUniform(S, kr, opt)
    S0 = kr.S;
    kr.S = S;
    kr.kr = interp1(S0, kr.kr, S);
end