function v = volumesOfTraps(Gt, res, traps, varargin)
% Compute volumes of (a subset of) precomputed structural traps
% 
% SYNOPSIS:
%   v = volumesOfTraps(Gt, res, varargin)
%
% PARAMETERS:
%   Gt       - top surface grid
%   res      - trap structure, as computed by trapAnalysis
%   traps    - vector with indices of traps for which to compute volumes.  If
%              empty, volumes are computed for _all_ traps.
%   varargin - keyword/value pairs for additional parameters.  Possible
%              options are:
%              * 'poro' - vector of porosities (default 1)
%              
%
% RETURNS:
%   v - vector containing volumes for each specified trap
%       NB: The volumes computes are total volumes, not pore volumes!
%
% SEE ALSO:
%   trapAnalysis    
%

    opt.poro = ones(Gt.cells.num, 1);
    opt = merge_options(opt, varargin{:});
    
    if isempty(traps)
       traps = 1:max(res.traps);
    end
    
    for i = 1:numel(traps)
        ind     = res.traps == traps(i);
        z       = Gt.cells.z(ind);
        fill_z  = res.trap_z(traps(i));

        assert(all((fill_z - z)>=0));

        H       = Gt.cells.H(ind);
        h_plume = min((fill_z - z),H);    
        v(i)    = sum(max(eps, Gt.cells.volumes(ind) .* h_plume .* opt.poro(ind)));
    end
end
