function [sW, sWcap] = getCapVisDist(block, pcVal, varargin)
% Upscale two-way oil-water distribution

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
    'gravity',  false  ... % include cravity in capillary part
	);
opt = merge_options(opt, varargin{:});

G  = block.G;
f  = block.fluid;
pv = block.pv;

if opt.gravity
    if isfield(f, 'rhoO')
        rhoO = f.rhoO;
    else
        rhoO = f.rhoOS;
    end
    if isfield(f, 'rhoW')
        rhoW = f.rhoW;
    else
        rhoW = f.rhoWS;
    end
    dRho = rhoW - rhoO;
    g    = 9.8066; % HARDCODED
    zi   = max(G.cells.centroids(:,3)) - G.cells.centroids(:,3);
    grav = -dRho.*g.*zi; % NOTE: Not sure about sign here

    sWcap = f.pcOWInv(pcVal - grav);
else
    % Saturation distribution from capillary limit
    sWcap = f.pcOWInv(pcVal.*ones(G.cells.num, 1));
end

% Divide domain into the x-direction
ijk  = gridLogicalIndices(G);
% [~,dir] = max(G.cartDims([1 2]));
% part = ijk{dir}; % 1->2 and 2->1
part = ijk{3};
numSubBlocks = numel(unique(part)); % Number of sub-blocks

% Use an upscaler to create sub blocks
upscaler = Upscaler(block.G, block.rock, 'deck', block.deck, ...
    'partition', part); % 'fluid', block.fluid, 

% Allocate space for distribution
sW = nan(G.cells.num,1);

% Loop over sub-blocks
for subInx = 1:numSubBlocks
    
    % Create grid, rock and fluid for sub block
    subCells = find(part==subInx);
        
    subBlock = upscaler.createBlock(subCells);
    subPv    = pv(subCells);
    subPvTot = sum(subPv);
    
    if subPvTot==0
        
        % The total pore volume is zero for the current sub block. This
        % must be handled separately. We simply leave the capillary
        % saturation values.
        subSw    = sWcap(subCells);
        
    else
        
        % Upscaled saturation for current sub block
        subSwU   = sum(sWcap(subCells).*subPv) / subPvTot;
        
        % Compute the saturation distribution using viscous limit 
        % upscaling, given the upscaled satuartion
        subSw    = getFracFlowDist(subBlock, subSwU);
        
    end
    
    % Store in vector
    sW(subCells) = subSw;
    
end

% TODO: TEMP Plotting for debugging purposes
if false
    cspan = [min(min(sW), min(sWcap)) max(min(sW), max(sWcap))];
    if cspan(1)==cspan(2)
        cspan = [0 1];
    end
    
    figure(1); clf;
    subplot(2,1,1); plotCellData(G, sWcap); view(0,0);
    colorbar; caxis(cspan);
    xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
    title('Capillary distribution');
    subplot(2,1,2); plotCellData(G, sW); view(0,0);
    colorbar; caxis(cspan);
    xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis');
    title('Capillary-viscous distribution');
    
    disp('hei');
end


end


function [sW, ff] = getFracFlowDist(block, sWU)
% Returnes the saturation distribution for a constant fractional flow value
% for the given upscaled water saturation.

pv    = block.pv;
pvTot = sum(pv);

% Create tables where T{i}=[ffW, sW].
T = createFracFlowTablesFromDeck(block.deck);

% Get all unique ff points for block
% As the relperm curves are piecewise linear functions, all the fractional 
% flow curves will also be piecewise linear functions between these points.
f = cellfun(@(x) x(:,1), T, 'UniformOutput', false);
f = unique([f{:}]);

% For each each fractional flow value, compute the corresponding upscaled
% water saturation. We compute all these points, because then we can do a
% simple 1D interpolation afterwards, and do not have to search for the
% solution.
reg    = handleRegions(block.deck);
satinx = reg.SATINX; % whole domain
if ischar(satinx)
    satinx = {satinx};
end
su = nan(size(f));
for i=1:numel(f)
    s     = interpReg(T, f(i).*ones(block.G.cells.num,1), satinx);
    su(i) = sum(s.*pv)/pvTot;
end

% Extend edges
% This ensures that interpolation below will not fail if sWU is numerically
% just outside of the interval su span
su = [su(1)-1; su; su(end)+1];
f  = [   f(1);  f;    f(end)];

% Compute the fractional flow value and the satruation distribution
% corresponding to the given upscaled saturation.
ff = interp1(su, f, sWU);
sW = interpReg(T, ff.*ones(block.G.cells.num,1), satinx);

% Check (as an implementation check)
diff = abs(sum(sW.*pv)/pvTot - sWU);
assert(diff < 1e-12, 'We are not close enough (diff=%1.2e)', diff);

end
