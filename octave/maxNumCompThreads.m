function res = maxNumCompThreads()
% Octave does not implement this MATLAB function.  This function provides and
% alternative.

res = nproc()/2; % return half the number of detected processors