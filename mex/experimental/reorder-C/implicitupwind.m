function state = implicitupwind(state, G, dt, rock, fluid, varargin)
%Compute solution of implicit upwind discretization twophase flow problem
%
% SYNOPSIS:
%   state = implicitupwind(state, G, dt, rock, fluid);
%   state = implicitupwind(state, G, dt, rock, fluid, 'src', src);
%   state = implicitupwind(state, G, dt, rock, fluid, 'src', src, ...
%                          'fluidopts', fluidopts, 'satnum', satnum);
%
%
% PARAMETERS:
%   state      - state struct with required fields 'flux' (vector of face
%                fluxes) and 's' (#cells x 2 array of phase saturations)
%
%   G          - grid structure.
%
%   dt         - time step or vector of sub time steps
%
%   rock       - rock struct with required field 'poro'.
%
%   dt         - vector of (sub)time steps to be computed
%
% OPTIONAL PARAMETERS
%   'src'      -
%
%   'bc'       -
%
%   'wells'    -
%
%   'satnum'   - double (or int32) vector of positive integers, one for each
%                cell, specifying which set of fluid parameters to use.
%                Numbering starts with one.
%
%   'fluidopt' - struct of parameters for corey fluid. Each field is a vector
%                with values for each saturation region:
%
%                'viscw' : water viscosity
%                'visco' : oil viscosity
%                'srw'   : residual water saturation
%                'sro'   : residual oil saturation
%                'nw'    : Corey exponent for water
%                'no'    : Corey exponent for oil
%
% RETURNS:
%   state      - Only the 's' field of state is set to the saturation for
%                each cell at the end of the time step(s).
%
% EXAMPLE:
%
% SEE ALSO:
%

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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


   % Initialize with 'params' field from fluid if possible
   if isfield(fluid, 'param'),
      opt  = struct('src', [], 'bc', [], 'wells', [], ...
                    'fluidopts', fluid.param,         ...
                    'satnum',    fluid.param.satnum,  ...
                    'substeps',  [],                  ...
                    'verbose',   false);
      opt.fluidopts = rmfield(opt.fluidopts, 'satnum');
   else
      opt  = struct('src', [], 'bc', [], 'wells', [], ...
                    'fluidopts', [], 'satnum', [],'substeps');
   end

   % Optional parameters take precedence
   opt  = merge_options(opt, varargin{:});

   check_fluidopts(opt);

   if ~isempty(opt.substeps),
      fprintf('\n\n');
      warning('''substeps'' option is deprecated, use vector of timesteps instead.');
      fprintf('\n\n');
   end

   flux = state.flux;
   s    = state.s(:,1);
   pv   = poreVolume(G, rock);
   q    = computeTransportSourceTerm(state, G, opt.wells, opt.src, opt.bc);
   q    = full(q);
   s    = implicitupwind_mex(G, pv, q, flux, dt, s, ...
                                int32(opt.satnum-1), opt.fluidopts);

   state.s = [s, 1-s];
end

function check_fluidopts(opt)
   if isempty(opt.fluidopts),
      error('No ''fluidopts'' found... unable to proceed');
   end

   if isempty(opt.satnum),
      error('No ''satnum'' found... using ones(...)');
      opt.satnum    = ones(G.cells.num, 1);
   end

   % Check existence of fileds in fluidopts
   assert(isfield(opt.fluidopts, 'viscw'), '''opt.fluidopts.viscw'' missing');
   assert(isfield(opt.fluidopts, 'visco'), '''opt.fluidopts.visco'' missing');
   assert(isfield(opt.fluidopts, 'srw'),   '''opt.fluidopts.srw'' missing');
   assert(isfield(opt.fluidopts, 'sro'),   '''opt.fluidopts.sro'' missing');
   assert(isfield(opt.fluidopts, 'nw'),    '''opt.fluidopts.nw'' missing');
   assert(isfield(opt.fluidopts, 'no'),    '''opt.fluidopts.no'' missing');

   % Check that max(satnum) match dimension of fields in fluidopts
   assert(numel(opt.fluidopts.viscw) == 1, ...
      'Number of entries in fluidopts.viscw must equal max(satnum).');
   assert(numel(opt.fluidopts.visco) == 1, ...
      'Number of entries in fluidopts.visco must equal max(satnum).');
   assert(numel(opt.fluidopts.srw) >= max(opt.satnum), ...
      'Number of entries in fluidopts.srw must equal max(satnum).');
   assert(numel(opt.fluidopts.sro) >= max(opt.satnum), ...
      'Number of entries in fluidopts.sro must equal max(satnum).');
   assert(numel(opt.fluidopts.nw) >= max(opt.satnum), ...
      'Number of entries in fluidopts.nw must equal max(satnum).');
   assert(numel(opt.fluidopts.no) >= max(opt.satnum), ...
      'Number of entries in fluidopts.no must equal max(satnum).');

end
