function [tau, pp] = timeStepping(time_param, print_param, iter)
% Calculates the next time step tau
%
% timeStepping determines the next time step according to the number of
% iterations from the last step in the following way:
%
%   If the number of iterations is less (or equal) than the lower optimal
%   range (lowerOptIterRange), it will increase the time step. In other words,
%   it will multiply the previous time step by a factor of lowerMultFactor.
%
%   On the other hand, if the number of iterations is greater (or equal)
%   than the upper optimal range (upperOptIterRange), it will decrease the
%   time step. Hence, it will multiply the previous time step by a
%   factor of upperMultFactor.
%                                                                                                                                                                                                                                                                            Finally, if the number of iterations lies within the lower and upper
%   If the time step lies within the optimal range, the previous
%   time step will remain unchanged.
%
% SYNOPSIS:
%   [tau, pp] = timeStepping(tau, tau_min, tau_max, simTime, ...
%                            timeCum, iter, printTimes, pp)
%
% PARAMETERS:
%   time_param   - Structure, containing time parameters
%   print_param  - Structure, containing print parameters
%   iter         - Scalar, number of iterations for the current time level
%
% RETURNS:
%   tau          - Scalar, time step for the next time level
%   pp           - Scalar, updated counter of printed times
%

%{
Copyright 2018-2020, University of Bergen.

This file is part of the fv-unsat module.

fv-unsat is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fv-unsat is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}

timeCum = time_param.time; % current time
tau_max = time_param.tau_max; % maximum time step
tau_min = time_param.tau_min; % minumum time step
simTime = time_param.simTime; % simulation time
dt_sim = time_param.tau; % auxiliar variable
lowerOptIterRange = 3; % lower optimal iteration range
upperOptIterRange = 7; % upper optimal iteration range
lowerMultFactor = 1.3; % lower multiplication factor
upperMultFactor = 0.7; % upper multiplication factor

printTimes = print_param.times; % printing times 
pp = print_param.print; % printing counter

% Time Step control
if (iter-1 <= lowerOptIterRange)
    dt_sim = dt_sim * lowerMultFactor;
    if (timeCum + dt_sim) > simTime %end of simulation
        tau = simTime - timeCum;
    elseif dt_sim > tau_max
        dt_sim = tau_max;
        tau = tau_max;
        if (timeCum + dt_sim) > printTimes(pp)
            tau = printTimes(pp) - timeCum;
            pp = pp +1;
            dt_sim = dt_sim / lowerMultFactor;
        end
    elseif (timeCum + dt_sim) > printTimes(pp)
        tau = printTimes(pp) - timeCum;
        pp = pp +1;
        dt_sim = dt_sim / lowerMultFactor;
    else
        tau = dt_sim;
    end
elseif (iter-1 >= upperOptIterRange)
    dt_sim = dt_sim * upperMultFactor;
    if (timeCum + dt_sim) > simTime %end of simulation
        tau = simTime - timeCum;
    elseif dt_sim < tau_min
        dt_sim = tau_min;
        tau = tau_min;
        if (timeCum + dt_sim) > printTimes(pp)
            tau = printTimes(pp) - timeCum;
            pp = pp +1;
            dt_sim = dt_sim * upperMultFactor;
        end
    elseif (timeCum + dt_sim) > printTimes(pp)
        tau = printTimes(pp) - timeCum;
        pp = pp +1;
        dt_sim = dt_sim * upperMultFactor;
    else
        tau = dt_sim;
    end
else
    if (timeCum + dt_sim) > simTime %end of simulation
        tau = simTime - timeCum;
    elseif (timeCum + dt_sim) > printTimes(pp)
        tau = printTimes(pp) - timeCum;
        pp = pp +1;
    else
        tau = dt_sim;
    end
end

end

