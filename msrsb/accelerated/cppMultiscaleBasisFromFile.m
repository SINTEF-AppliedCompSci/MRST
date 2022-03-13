function I = cppMultiscaleBasisFromFile(fn)
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
    [Nf, Nc, fine, coarse] = readGrid(fn);
    
    out = fullfile(fn, 'output');
    if ~exist(out, 'dir')
        mkdir(out)
    end
    
    mex_cppMultiscaleBasis(fn);
    
    I_comp = readOperator(out);
    
    I = sparse(fine, coarse, I_comp, Nf, Nc);
end


function [Nf, Nc, fine, coarse] = readGrid(fn)
    inp = fullfile(fn, 'input');
    info = dlmread(fullfile(inp, 'info.txt'));
    
    Nf = info(1, 1);
    Nc = info(2, 1);
    offsets = info(3, :) + 1;
    fine = dlmread(fullfile(inp, 'support.txt')) + 1;
    coarse = rldecode((1:Nc), diff(offsets), 2);
end

function I_comp = readOperator(out)
    I_comp = dlmread(fullfile(out, 'operator.txt'));
end