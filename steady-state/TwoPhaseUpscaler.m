classdef TwoPhaseUpscaler < OnePhaseUpscaler
    %Two phase upscaling

properties
    RelpermMethod     % Relperm upscaling method
    RelpermAbsMethod  % Abs-perm upscaling used in relperm upscaling
    nrelperm
    pcow
    npcow
    pcowgrav     % Whether to include gravity in pcOW upscaling
    
    relpermdims  % Dimensions to upscale relperm in
    savesat      % Save saturation distributions
    
end

methods
    
    function upscaler = TwoPhaseUpscaler(G, rock, fluid, varargin)
        upscaler = upscaler@OnePhaseUpscaler(G, rock, 'fluid', fluid, ...
            varargin{:});
        
        upscaler.RelpermMethod = [];
        upscaler.RelpermAbsMethod = 'pressure';
        upscaler.nrelperm = 20;
        upscaler.pcow     = true; % Upscale pcow or not
        upscaler.npcow    = 100;
        upscaler.pcowgrav = true;
        upscaler.relpermdims = upscaler.dims;
        upscaler.savesat  = false; % Save saturation distributions
        
        upscaler = merge_options(upscaler, varargin{:});
    end
    
    function [blockdata, globdata, report] = upscale(upscaler)
        
        [blockdata, globdata, report] = upscale@OnePhaseUpscaler(upscaler);
        
        % Special treatment of the endpoint scaling of two-phase
        if strcmpi(upscaler.RelpermMethod, 'endpoint')
            % Once we are done with the upscaling of each block, we perform
            % an average of the relperm and pcow
            fprintf(['Computing global average of two-phase '...
                'properties...\n']);
            b = GridBlock(upscaler.G, upscaler.rock, ...
                'deck', upscaler.deck, 'fluid', upscaler.fluid);
            [d, r] = upRelPermPV(b, [], 'porevolume', ...
                'nvalues', upscaler.nrelperm);
            globdata.relperm      = d;
            report.global.relperm = r;
        end
        
    end
    
    function [data, report] = upscaleBlock(upscaler, block)
        
        wantReport  = nargout > 1;
        startTime   = tic;
        
        % Perform one phase upscaling first
        [data, rep] = upscaleBlock@OnePhaseUpscaler(upscaler, block);
        if wantReport
            report.onephase = rep;
        end
        
        data.relpermdims = upscaler.relpermdims;
        
        % Capillary pressure upscaling
        % If using pore volume relperm upscaling, we also upscale the pcOW
        % in the relperm function, so don't do it here.
        if upscaler.pcow && ~any(strcmpi(upscaler.RelpermMethod, ...
                {'porevolume', 'endpoint'}))
            
            if ~isempty(regexpi(upscaler.RelpermMethod, '_grav$'));
                
                if upscaler.pcowgrav
                    % We upscale with gravity
                    [data, rep] = upPcOW(block, data, ...
                        'npointsmax', upscaler.npcow, ...
                        'gravity', 'centroid');
                else
                    % We upscale without gravity
                    [data, rep] = upPcOW(block, data, ...
                        'npointsmax', upscaler.npcow);
                end
                
                % We also need the bottom pcow curve for relperm upscaling
                [d, r] = upPcOW(block, [], ...
                    'npointsmax', upscaler.npcow, 'gravity', 'bottom');
                data.pcOW_bot = d.pcOW;
                if wantReport
                    report.pcow_bot = r;
                end
            else
                % Default: no gravity
                [data, rep] = upPcOW(block, data, ...
                    'npointsmax', upscaler.npcow);
            end
            if wantReport
                report.pcow = rep;
            end
            if upscaler.verbose
                fprintf('  Cap.pres:     %6.3fs\n', rep.time);
            end
        end
        
        % Relative permeability upscaling
        if strcmpi(upscaler.RelpermMethod, 'porevolume')
            up = @() upRelPermPV(block, data, upscaler.RelpermMethod, ...
                'nvalues', upscaler.nrelperm);
        elseif strcmpi(upscaler.RelpermMethod, 'endpoint')
            % Enpoint scaling of SWOF data
            up = @() upRelPermEPS(block, data, upscaler.RelpermMethod, ...
                'absmethod', upscaler.RelpermAbsMethod, ...
                'dp', upscaler.dp, 'dims', upscaler.relpermdims);
        elseif strcmpi(upscaler.RelpermMethod, 'endpointfull')
            % Enpoint scaling of SWOF data, but with SWOF upscaling in each
            % grid block
            up = @() upRelPermEPS(block, data, upscaler.RelpermMethod, ...
                'absmethod', upscaler.RelpermAbsMethod, ...
                'dp', upscaler.dp, 'dims', upscaler.relpermdims, ...
                'fullswof', true);
        else
            up = @() upRelPerm(block, data, upscaler.RelpermMethod, ...
                'nsat', upscaler.nrelperm, 'dp', upscaler.dp, ...
                'absmethod', upscaler.RelpermAbsMethod, ...
                'dims', upscaler.relpermdims, 'savesat', upscaler.savesat);
        end
        
        % Upscale
        if wantReport
            [data, rep] = up();
            report.relperm = rep;
        else
            data = up();
        end
        
        if upscaler.verbose
            fprintf('  Rel.perm:     %6.3fs\n', rep.time);
        end
        
        if wantReport
            totalTime   = toc(startTime);
            report.time = totalTime;
        end
        
    end
    
    function block = createBlock(upscaler, cells)
        % Create grid, rock and fluid for the sub block represented by the
        % given cells.
        
        block = createBlock@OnePhaseUpscaler(upscaler, cells);
        
        if ~isempty(block.deck)
            
            % To use capillary upscaling, we need a function pcOWInv
            block.fluid = addPcOWInvADIFluid(block.fluid, block.deck);
            
            % To use viscous upscaling, we need a function fracFlowInv
            block.fluid = addFracFlowInvADIFluid(block.fluid, block.deck);
            
        end
        
    end
    
end

end

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
