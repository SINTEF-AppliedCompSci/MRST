function varargout = msfvCPRWrapper(fn, varargin)
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

    persistent itercount oldsat;
    itercount = 0;

    varargout = cell(nargout, 1);
    [varargout{:}] = fn(varargin{:});


    state = varargout{1};
    if isempty(oldsat);
        oldsat = state.s;
    end

    if varargout{2}.converged
        itercount = itercount + 1;
%         updateTargets = sum(abs(oldsat - state.s) > 0.2, 2) > 0;
%         if any(updateTargets)
%             msfvPressureSolve(updateTargets);
%         end


%         if itercount == 1
          if state.T > 10*day

            % Reset solver
            disp('Converged, resetting fv basis!')
            msfvPressureSolve();
            itercount = 0;
            state.T = 0;
        end
    end
end
