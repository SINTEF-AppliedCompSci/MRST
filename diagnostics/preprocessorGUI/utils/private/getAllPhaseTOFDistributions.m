function [data,tof,tof_ix] = getAllPhaseTOFDistributions(state, pv, wellIx, wellType, D, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('maxTOF',     500*year,...
                 'normalize',  'true');
    opt = merge_options(opt, varargin{:});
    if isempty(opt.maxTOF)
        opt.maxTOF = inf;
    end
    
    switch wellType
        case 'producer'
            tofNo = 2;
            iReg  = D.ptracer;
        case 'injector'
            tofNo = 1;
            iReg  = D.itracer;
        otherwise
            error('Incorrect well type');
    end
    
    tof = D.tof(:,tofNo)./year;    
    % Don't include maxTOF (avoid jump at maxTOF in plot)
    tof_ix = find(and(iReg(:,wellIx)>.01,  tof < .99*opt.maxTOF/year ));
    
    % Sort tof
    [tof, ix] = sort(tof(tof_ix));
    tof_ix = tof_ix(ix);
    
    % Compute fluid mobilities
    if isfield(state, 'mob')
        % Use fraction of total mobility to estimate how much will flow
        data = bsxfun(@rdivide, state.mob, sum(state.mob, 2));
    else
        data = state.s;
    end
    
    % Weight by pore volumes*tracervalue to get actual volumes 
    data = bsxfun(@times, data, pv.*iReg(:,wellIx));
    data = data(tof_ix,:);
    
    % Cumsum data
    data = cumsum(data);
    
    % use approx 50 points in plot
    di   = ceil(numel(tof)/50);
    tof    = tof(1:di:end);
    tof_ix    = tof_ix(1:di:end);
    data = data(1:di:end, :);
    if opt.normalize
       data = bsxfun(@rdivide, data, sum(data, 2));
    end
    data(isnan(data)) = 0;
end

