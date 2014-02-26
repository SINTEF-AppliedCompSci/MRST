function v = volumesOfTraps(Gt, res, varargin)
%% Compute volumes of (a subset of) precomputed structural traps
% NB: The volumes computes are total volumes, not pore volumes! 
% 
% SYNOPSIS:
%   function v = volumesOfTraps(Gt, res, varargin)
%
% PARAMETERS:
%   Gt       - top surface grid
%   res      - trap structure, as computed by trapAnalysis
%   varargin - may contain a vector with indices of traps for which to
%              compute volumes, or can be left empty, if volumes are to be
%              computed for _all_ traps
%
% RETURNS:
%   v - vector containing volumes for each specified trap
%
% SEE ALSO:
% trapAnalysis    
%

    if nargin == 2
        trap = 1:max(res.traps);
    else
        trap = varargin{1};
    end
    v = zeros(1,numel(trap));
    for i = 1:numel(trap)
        ind     = res.traps == trap(i);
        z       = Gt.cells.z(ind);
        fill_z  = res.trap_z(trap(i));

        assert(all((fill_z - z)>=0));

        H       = Gt.cells.H(ind);
        h_plume = min((fill_z - z),H);    
        v(i)    = sum(max(eps, Gt.cells.volumes(ind).*h_plume));
    end
end
