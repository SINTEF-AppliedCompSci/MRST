function dispControls(controls, schedule, varargin)
%Display control values.

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


fprintf('\n----------------- DISPLAYING CONTROL VARIABLES ----------------\n')

numC    = numel(controls.well);
cw      = [controls.well(:).wellNum];
names   = [schedule(1).names(cw)];
types   = {controls.well(:).type};
minMax  = vertcat(controls.well.minMax);

fprintf('%9s%9s%9s%15s\n', 'Var', 'Name', 'Type', 'MaxMin')
for k = 1: numC
    var     = sprintf('u_%d', k);
    vars{k} = var;
    fprintf('%9s%9s%9s%15s\n', var, names{k}, types{k}, mat2str(minMax(k, :)) );
end

ec = controls.linEqConst;
if ~isempty(ec)
    fprintf('\nLinear equality constraints: \n')
    numEC = size(ec.A, 1);
    for k = 1: numEC
        A = ec.A(k, :); b = ec.b(k);
        [r, c] = size(A);
        lstart = false;
        for k1 = 1:c
            if A(k1) ~= 0
                if lstart, fprintf(' + '); end
                if A(k1) ~= 1, fprintf(num2str(A(k1))); end
                fprintf(vars{k1});
                lstart = true;
            end
        end
        fprintf(' = ');
        fprintf(num2str(b));
        fprintf('\n');
    end
end




