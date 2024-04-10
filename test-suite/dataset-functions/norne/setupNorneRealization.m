function  [G, rock, fluid, deck] = setupNorneRealization(realNo, varargin)
% Helper function to generate specific Norne realization. Will run 
% generateNorneEnsemble.m if this hasn't been done already. 
%
% SYNOPSIS:
%   out = setupNorneRealization(realNo)
%
% DESCRIPTION:
%    This function returns the variables needed to setup a model of a
%    specific Norne realization.
%
% PARAMETERS:
%    realNo  - requested realization. Number between 1 and 50. If omitted,
%              defaults to 1.
%
% KEYWORD ARGUMENTS:
%    actnum  - override for ACTNUM. Must be a subset of the ACTNUM in the
%              generated ensemble. This is checked by the function.
%
% RETURNS:
%     G - MRST grid structure for realization realNo.
%     rock - rock structure for realization realNo.
%     fluid - MRST fluid structure.
%     deck - deck structure for realization realNo.
%    
% EXAMPLE:
%   [G, rock, fluid, deck] = setupNorneRealization(realNo)
%
% SEE ALSO:
%   `generateNorneEnsemble`, `setupNorneFn`

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

    if nargin<1    
        realNo = 1;
    end
    opt = struct('actnum', []);
    opt = merge_options(opt, varargin{:});
  
    try
        pth = getDatasetPath('norne_ensemble');
        norne=load(fullfile(pth,'data','NorneInitEns.mat'));
    catch me
        if ~exist('pth', 'var')
           error('Download:Failed', ...
                 'Failed to download Norne Ensemble dataset: %s', ...
                 me.message);
        end
        addpath(fullfile(pth,'data'));
        generateNorneEnsemble(max(50,realNo));
        rmpath(fullfile(pth,'data'));
        
        norne = load(fullfile(pth,'data','NorneInitEns.mat'));
    end
    assert(realNo<=numel(norne.ensemble));
    data = norne.ensemble(realNo);

    % The generated data contains 12 extra disconnected cells that need to
    % be accounted for when assigning data
    if ~ (makeNorneSubsetAvailable() && makeNorneGRDECL())
        error('Unable to obtain simulation model subset');
    end
    
    grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');
    grdecl = convertInputUnits(readGRDECL(grdecl), getUnitSystem('METRIC'));
    actnum = opt.actnum;
    if ~isempty(actnum)
        assert(numel(actnum) == numel(grdecl.ACTNUM));
        imap0 = find(grdecl.ACTNUM);
        imap = find(actnum);
        overlap = setdiff(imap, imap0);
        assert(isempty(overlap));
        % TODO: Is ordering ok here?
        newix = ismember(imap, imap0);
        data.PERMX = data.PERMX(newix);
        data.PORO = data.PORO(newix);
        data.NTG = data.NTG(newix);
        grdecl.ACTNUM = actnum;
    end
    G = processGRDECL(grdecl, 'SplitDisconnected',false);
    
    % Assign permeability: apply multipliers as in ECLIPSE input file
    data.PERMX = convertFrom(exp(data.PERMX), milli*darcy);
    data.PERMY = data.PERMX;

    data.PERMZ = data.PERMX.*grdecl.PERMZ(G.cells.indexMap)./grdecl.PERMX(G.cells.indexMap);

    % Assign additional multipliers from the ensemble data
    data.MULTZ = grdecl.MULTZ(G.cells.indexMap);
    [~,~,K] = gridLogicalIndices(G);
    n=0;
    for k=[1,8,11,12,15,18]
        ind = (K==k);
        nk = sum(ind);
        data.MULTZ(ind) = 10.^data.multz(n+1:n+nk);
        n = n + nk;
    end
    
    % Reconstruct the grid allowing the disconnected parts to be split
    g      = processGRDECL(grdecl);
    G      = computeGeometry(g(1));

    % Disregard g(2) when extracting the petrophysical data
    imap   = [g(1).cells.indexMap; g(2).cells.indexMap];
    actnum = [true(g(1).cells.num,1); false(g(2).cells.num,1)];
    [~,i]  = sort(imap);

    ind  = find(actnum(i));
    rock = grdecl2Rock(data, ind);
    deck.GRID.MULTZ = data.MULTZ(ind);

    % Set fluid data
    fluid = initSimpleADIFluid('phases', 'WO', 'mu', [1,1]*centi*poise, 'n', [2,2]);
end
