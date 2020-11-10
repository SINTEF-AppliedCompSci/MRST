function simulateEnsembleMembersStandalone(fileName, range)
% Utility function for running ensemble simulations in a background
% process.
%
% SYNOPSIS:
%   simulateEnsembleMembersStandalone(fileName, range)
%
% PARAMETERS:
%   fileName - a .mat file with a stored 'ensemble' object.
%   range    - the ensemble members that should be simulated.

%{
#COPYRIGHT#
%}
    data = load(fileName);
    data.ensemble.simulateEnsembleMembersCore(range);
end