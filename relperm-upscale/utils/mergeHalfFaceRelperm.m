function kr = mergeHalfFaceRelperm(model, kr, varargin)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

    opt = struct('type', 'cell');
    [opt, extra] = merge_options(opt, varargin{:});
    
    for ph = 1:numel(kr)
        switch lower(opt.type)
            case 'cell'
                n = model.G.cells.num;
                mapping = model.operators.N;
                kr_next = cell(n, 1);
            case 'face'
                n = size(model.operators.N, 1);
                mapping = repmat((1:n)', 1, 2);
                kr_next = cell(n, 1);
            otherwise
                error('Unknown merge strategy');
        end
        kr_old = kr{ph}.reservoir;
        mapping = mapping(:);
        for i = 1:n
            ix = (mapping == i);
            
            kr_next{i} = mergeTables(kr_old{ix});
        end        
        kr{ph}.reservoir = kr_next;
    end
end

function t = mergeTables(varargin)
    S = [];
    kr = [];
    for i = 1:nargin
        S = [S; varargin{i}.S];
        kr = [kr; varargin{i}.kr];
    end
    [S, ix] = sort(S);
    kr = kr(ix);
    
    t.S = S;
    t.kr = kr;
end
