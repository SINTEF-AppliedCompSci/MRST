function rock = SPE10_rock(varargin)
%Define rock properties for Model 2 of tenth SPE CSP
%
% SYNOPSIS:
%   rock = SPE10_rock
%   rock = SPE10_rock(layers)
%   rock = SPE10_rock(I, J, K)
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
%   The permeability data is returned unconverted, meaning the data is
%   given in whatever units are provided in the MAT-file representation
%   "spe10_rock.mat" in the directory containing "SPE10_rock.m".
%   This is typically milli*darcy.
%
%   It is the caller's responsibility to convert this data into MRST's
%   strict SI-only unit conventions.  Function 'convertFrom' provides one
%   conversion method.
%
% EXAMPLE:
%   rock = SPE10_rock(85)
%
% SEE ALSO:
%   SPE10_setup, make_spe10_data, convertFrom, milli, darcy.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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

   [Nx, Ny, Nz] = deal(60, 220, 85);
                                          % ()
   [I, J, K]    = deal(1:Nx, 1:Ny, 1:Nz); % Default (entire dataset)

   if nargin == 1,
                                          % (layers)
      K = varargin{1};                    % Caller specified ind. layers

   elseif nargin == 3,
                                          % (I, J, K)
      [I, J, K] = deal(varargin{:});      % Caller specified box

   elseif nargin ~= 0,
      error(['Syntax is\n\t'              , ...
             'rock = %s         %% or\n\t', ...
             'rock = %s(layers) %% or\n\t', ...
             'rock = %s(I, J, K)'], mfilename, mfilename, mfilename);
   end

   assert (all((0 < I) & (I <= Nx)), 'I outside valid range 1:60');
   assert (all((0 < J) & (J <= Ny)), 'J outside valid range 1:220');
   assert (all((0 < K) & (K <= Nz)), 'K outside valid range 1:85');

   % Load data in pre-processed form.
   rock_file = fullfile(fileparts(mfilename('fullpath')), 'spe10_rock');
   if ~exist([rock_file, '.mat'], 'file'),
      ok = make_spe10_data;
      assert (ok);
   end
   data = load(rock_file);
   rock = data.rock; % Fa-Fa-Fa

   % Extract caller's requested subset from 'rock' data
   %
   % Return only PERM and PORO data and exclude any other information that
   % might be stored in on-disk representation of rock data.
   %
   [I, J, K] = ndgrid(I, J, K);
   ix = sub2ind([Nx, Ny, Nz], I(:), J(:), K(:));

   rock = struct('perm', rock.perm(ix, :), ...
                 'poro', rock.poro(ix   ));
end
