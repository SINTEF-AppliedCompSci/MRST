function rock = getSPE10rock(varargin)
%Define rock properties for Model 2 of tenth SPE CSP
%
% SYNOPSIS:
%   rock = getSPE10rock
%   rock = getSPE10rock(layers)
%   rock = getSPE10rock(I, J, K)
%
% PARAMETERS:
%   layers - Which of the 85 model layers to include in a specific test.
%            OPTIONAL.  If unspecified, all 85 layers (for a total of
%            60-by-220-by-85 == 1,122,000 cells) are included.
%
%            Some possible values are
%               layers = ( 1 : 35).';  %  Tarbert formation
%               layers = (36 : 85).';  %  Upper Ness formation
%
%   I,J,K  - Global Cartesian indices identifying subset of rock data.
%            Useful, e.g., in extracting lateral sections.  Arrays must
%            satisfy
%
%               ALL((0 < I) & (I <=  60)) && ...
%               ALL((0 < J) & (J <= 220)) && ...
%               ALL((0 < K) & (K <=  85))
%
%            OPTIONAL.  If unspecified, treated as I=1:60,J=1:220,K=1:85
%            (i.e., the entire model; 1,122,000 cells).  In other words,
%            equivalently to an unspecified 'layers' input.
%
% RETURNS:
%   rock - Rock structure having fields 'perm' and 'poro' pertaining to
%          the specified layers or cell subset.
%
% NOTE:
%   The permeability data is returned in strict SI units (metres squared).
%
% EXAMPLE:
%   rock = getSPE10rock(85)
%
% SEE ALSO:
%   `getSPE10setup`, `make_spe10_data`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

   assert (all(cellfun(@isnumeric, varargin)), ...
           'Input arrays must be numeric.');

   % Load data in pre-processed form.
   rock = load_mat_file();

   % Extract caller's requested subset from 'rock' data
   %
   % Return only PERM and PORO data and exclude any other information that
   % might be stored in on-disk representation of rock data.
   %
   ix   = define_subset(varargin{:});

   rock = struct('perm', rock.perm(ix, :), ...
                 'poro', rock.poro(ix   ));
end

%--------------------------------------------------------------------------

function ix = define_subset(varargin)
   [Nx, Ny, Nz] = deal(60, 220, 85);
                                          % ()
   [I, J, K]    = deal(1:Nx, 1:Ny, 1:Nz); % Default (entire dataset)

   if nargin == 1
                                          % (layers)
      K = varargin{1};                    % Caller specified ind. layers

   elseif nargin == 3
                                          % (I, J, K)
      [I, J, K] = deal(varargin{:});      % Caller specified box

   elseif nargin ~= 0
      file = mfilename();
      error(['Syntax is\n\t'              , ...
             'rock = %s         %% or\n\t', ...
             'rock = %s(layers) %% or\n\t', ...
             'rock = %s(I, J, K)'], file, file, file);
   end

   validate_range(I, Nx);
   validate_range(J, Ny);
   validate_range(K, Nz);

   [I, J, K] = ndgrid(I, J, K);

   ix = sub2ind([Nx, Ny, Nz], I(:), J(:), K(:));
end

%--------------------------------------------------------------------------

function rock = load_mat_file()
   rdir = getDatasetPath('spe10', 'skipAvailableCheck', true);
   rock_file = fullfile(rdir, 'spe10_rock');

   if ~exist([rock_file, '.mat'], 'file')
      ok = makeSPE10DataAvailable();

      if ~ok
         error('SPE10Download:Fail', 'Failed to download SPE 10 dataset');
      end
   end

   data = load(rock_file);
   rock = data.rock; % Fa-Fa-Fa
end

%--------------------------------------------------------------------------

function validate_range(i, n)
   assert (all((0 < i) & (i <= n)), ...
           '%s outside valid range 1:%d', inputname(1), n);
end
