function schedule = simple_injection_migration_schedule(W, bc, inj_duration, ...
                                                        inj_steps, migration_duration, ...
                                                        migration_steps)

    % Define simple, single-well injection and migration scenario

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

    % Setting up two copies of the well and boundary specifications. 
    % Modifying the well in the second copy to have a zero flow rate.
    schedule.control    = struct('W', W, 'bc', bc);
    schedule.control(2) = struct('W', W, 'bc', bc);
    schedule.control(2).W.val = 0;

    dT_injection = rampupTimesteps(inj_duration, ...
                                   inj_duration/inj_steps, 7);
   
    dT_migration = repmat(migration_duration/migration_steps, migration_steps, 1);
    
    schedule.step.val = [dT_injection; dT_migration];
    schedule.step.control = [ones(numel(dT_injection), 1); ...
                             2 * ones(numel(dT_migration), 1)];

end
