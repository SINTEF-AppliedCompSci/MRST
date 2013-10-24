function command = simcommand(dir, sim)

% Construct a simulator command depending on user and context.
%
% SYNOPSIS:
%   command = simcommand(context, sim, user)
%
% PARAMETERS:
%   context - can be either 'opm' or 'urc'
%   sim     - the filename of the simulator you want to run,
%             including necessary paths from the root, for
%             example 'examples/sim_2p_comp_reorder'
%   user    - will be passed to opm_dir() or urc_dir()
%             depending on context
%
% RETURNS:
%   A string with the full path to the requested simulator program.

%{
#COPYRIGHT#
%}

% This seems necessary on some Linux variants, and does no harm on Mac OS.
ldfix='export LD_LIBRARY_PATH=/lib/x86_64-linux-gnu/:/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH; ';

command = [ldfix, fullfile(dir, sim)];
