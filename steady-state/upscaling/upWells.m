function CW = upWells(CG, rock, W, varargin)
% Upscale wells

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

opt = struct(...
    'LinSolve',  @mldivide, ...
    'bc_method',  'wells_simple', ...
    'debug',     false ...
    );
opt = merge_options(opt, varargin{:});

if opt.debug
    % We are in debug mode. We do not spend time on performing well
    % upscaling. Instead, we just add some dummy well indecies.
    CW = makeDummyCGWells(CG, W);
    return;
end

% Compute transmissibility
s = setupSimComp(CG.parent, rock);
T = s.T_all;

% Perform upscaling
warning('off','upscaleTransNew:ZeroBoundary');
[~, ~, CW, ~] = upscaleTransNew(CG, T, ...
   'wells', {W}, 'bc_method', opt.bc_method, ...
   'LinSolve', opt.LinSolve );
warning('on','upscaleTransNew:ZeroBoundary');

end


% Function only used for debugging purposes. Creates a dummy upscaled well
% structure with dummy upscaled well indecies. The method is copied from
% upscaleTransNew, and just slightly altered.
function cgwells = makeDummyCGWells(cg, wells)
cgwells = wells;
for i = 1 : numel(wells),
    fcells  = wells(i).cells;
    cgcells = cg.partition(fcells);
    
    tab        = sortrows([cgcells, fcells, (1 : numel(fcells)) .']);
    [cells, n] = rlencode(tab(:,1));
    fcellspos  = cumsum([1 ; n]);
    
    if cg.griddim > 2,
        pno = rldecode(1 : numel(cells), n, 2) .';
        cc  = cg.parent.cells.centroids(tab(:,2), 3);
        cv  = cg.parent.cells.volumes  (tab(:,2));
        
        hpos = sparse(pno, 1 : numel(pno), cv) ...
            * [ cc, ones([numel(pno), 1]) ];
        
        hpos = hpos(:,1) ./ hpos(:,2);         clear pno cc cv
    else
        hpos = 0;
    end
    
    % Compute WI as some simple average. This is instead of actually
    % upscaling the WI value, and is only used for debugging purposes.
    wi = sum(wells(i).WI.*wells(i).dZ)./sum(wells(i).dZ);
    
    cgwells(i).cells     = cells;
    cgwells(i).WI        = wi.*ones([numel(cells), 1]);
    cgwells(i).dZ        = hpos - wells(i).refDepth;
    cgwells(i).fcellspos = fcellspos;
    cgwells(i).fcells    = tab(:,2);
    cgwells(i).fperf     = tab(:,3);
end

cgwells(i).parent = wells(i);
end

