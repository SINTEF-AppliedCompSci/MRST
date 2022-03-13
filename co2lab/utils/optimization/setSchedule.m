function schedule = setSchedule(Gt, rock, wcells, qtot, isteps, itime, ...
                                    msteps,mtime, single_control, varargin) 
% Construct a schedule object that is convenient for use with
% 'optimizeFormation'.  Wells are rate-controlled.
%
% SYNOPSIS: 
% function schedule = setSchedule(Gt, rock, wcells, qtot, isteps, itime, ...
%                                 msteps,mtime, single_control, varargin)
%
% PARAMETERS:
%   Gt             - the simulation grid (top surface grid)
%   rock           - associated rock structure
%   wcells         - well cells (vector containing the indices of the well
%                    cells, one entry per welll)
%   qtot           - total amount to inject into each well (vector with one
%                    entry per well).  Rates will be obtained by dividing
%                    entries of this vector by total injection time.
%   isteps         - number of injection time steps 
%   itime          - duration of injection phase
%   msteps         - number of migration time steps
%   mtime          - duration of migration phase
%   single_control - whether the injection phase should be governed by one
%                    single control, or a new control for each timestep
%                    (allowing for dynamic adjustment of rates)
%   varargin       - The only optional argument is 'minval', which specifies
%                    the minimum possible rate for a well.
%
% RETURNS:
%   schedule - the constructed schedule
%
% SEE ALSO:
%  `optimizeFormation`

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

   opt.minval = 0;
   opt = merge_options(opt, varargin{:});
   
    assert(isteps>0);
    if msteps == 1
        msteps = 2;
        warning(['If migration is happening, we need at least two steps, in ' ...
                 'order to have a transition step.  Migration step has been ' ...
                 'increased to two.']);
    end
    
    % constructing wells
    W = [];
    wellradius = 0.3;
    
    % computing fixed rates
    rate = qtot / itime;  
    
    for i = 1:numel(wcells)
       W = addWellVE(W, Gt, rock, wcells(i) , ...
                     'Type'   , 'rate'      , ...
                     'Val'    , rate(i)     , ...
                     'Radius' ,  wellradius , ...
                     'Comp_i' , [0, 1]      , ...
                     'name'   , ['I', num2str(i)]);
    end

    % constructing schedule
    cpos = 1;
    ctrls = [];%#ok
    if single_control && isteps > 0 
        schedule.control(cpos).W = W;
        cpos = cpos+1;
        ctrls = ones(isteps,1);
    else
        for i = 1:isteps
            schedule.control(cpos).W = W;
            cpos = cpos+1;
        end
        ctrls = [1:isteps]';%#ok
    end    
    
    if msteps > 0
        schedule.control(cpos).W = W;
        for i = 1:numel(schedule.control(cpos).W)
            schedule.control(cpos).W(i).val = opt.minval;
        end
        cpos = cpos+1;%#ok
    end

    mig_ctrl = numel(schedule.control); % control step for migration, only
                                        % relevant (and correct) if msteps > 0

    dTi = itime/isteps;
    dTm = mtime/(msteps-1);
    
    istepvec = ones(isteps, 1) * dTi; 
    mstepvec = [dTi; ones(msteps-1,1) * dTm];
    mstepvec(2) = mstepvec(2) - dTi;

    schedule.step.val = [istepvec; mstepvec];
    schedule.step.control = [ctrls; mig_ctrl * ones(msteps,1)];
end
