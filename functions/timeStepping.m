function [dt,printCounter] = timeStepping(dt,dt_min,dt_max,simTime,timeCum,iter,printTimes,printCounter,...
                                 lowerOptIterRange,upperOptIterRange,lowerMultFactor,upperMultFactor)
% Calculates the next time step dt
%
% timeStepping determines the next time step according to the number of 
% iterations from the last step in the following way.
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
%   If the time step lies withing the optimal range, the previous 
%   time step will remain unchanged.
%
% SYNOPSIS:
%   [dt,printCounter] = timeStepping(dt,dt_min,dt_max,simTime,timeCum,iter,printTimes,printCounter,...
%                                    lowerOptIterRange,upperOptIterRange,lowerMultFactor,upperMultFactor)
%
% PARAMETERS:
%   dt                  - Scalar, previous time step
%   dt_min              - Scalar, minimum time step                 
%   dt_max              - Scalar, maximum time step
%   simTime             - Scalar, final simulation time
%   timeCum             - Scalar, current simulation time
%   iter                - Scalar, number of iterations of the current time
%   printTimes          - Vector, containing printing times
%   printCounter        - Scalar, counter of printed times
%   lowerOptIterRange   - Scalar, lower optimal iteration range
%   upperOptIterRange   - Scalar, upper optimal iteration range
%   lowerMultFactor     - Scalar, lower multiplication factor
%   upperMultFactor     - Scalar, upper multiplication factor
%
%  RETURNS:
%   dt                  - Scalar, time step for the next time level
%   printCounter        - Scalar, updated counter of printed times
%

dt_sim = dt;
c_low = lowerMultFactor;
c_upp = upperMultFactor;

% Time Step control
    if (iter-1 <= lowerOptIterRange)
        dt_sim = dt_sim * c_low;
        if (timeCum + dt_sim) > simTime %end of simulation
            dt = simTime - timeCum;
        elseif dt_sim > dt_max
            dt_sim = dt_max;
            dt = dt_max;
            if (timeCum + dt_sim) > printTimes(printCounter)
                dt = printTimes(printCounter) - timeCum;
                printCounter = printCounter +1;
                dt_sim = dt_sim / c_low;
            end
        elseif (timeCum + dt_sim) > printTimes(printCounter)
            dt = printTimes(printCounter) - timeCum;
            printCounter = printCounter +1;
            dt_sim = dt_sim / c_low;
        else
            dt = dt_sim;
        end
    elseif (iter-1 >= upperOptIterRange)
        dt_sim = dt_sim * c_upp;
        if (timeCum + dt_sim) > simTime %end of simulation
            dt = simTime - timeCum;
        elseif dt_sim < dt_min
            dt_sim = dt_min;
            dt = dt_min;
            if (timeCum + dt_sim) > printTimes(printCounter)
                dt = printTimes(printCounter) - timeCum;
                printCounter = printCounter +1;
                dt_sim = dt_sim * c_upp;
            end
        elseif (timeCum + dt_sim) > printTimes(printCounter)
            dt = printTimes(printCounter) - timeCum;
            printCounter = printCounter +1;
            dt_sim = dt_sim * c_upp;
        else
            dt = dt_sim;
        end
    else
        if (timeCum + dt_sim) > simTime %end of simulation
            dt = simTime - timeCum;
        elseif (timeCum + dt_sim) > printTimes(printCounter)
            dt = printTimes(printCounter) - timeCum;
            printCounter = printCounter +1;
        else
            dt = dt_sim;
        end
    end
end

