function deck = refineDeck(deck_in, dim, varargin)
%
% Refine the grid resolution of a deck, and update other information
% (`RUNSPEC`, wells, selected cell-based fields under `SOLUTION`, `REGIONS`
% and `PROPS`) accordingly. 
% 
% Cell-based fields are updated by letting refined cells inherit the values from
% their 'parent' coarse cells.  Wells are updated by changing the cell indexing,
% adding perforations in the direction of the original perforations, and divide
% refined well indices and KH by the refinement in the perforation direction.
%
% NOTE:
%   This function is not fully tested and has only been used on a limited
%   subset of models. While potentially useful, it should be used with care
%   and results should be carefully examined.
%
% SYNOPSIS:
%   function deck = refineDeck(deck_in, dim, varargin)
%
% PARAMETERS:
%   deck_in  - deck to refine
%   dim      - 3-component vector with refinement factor in each cooridinate
%              direction
%   varargin - supports the (true/false) keyword 'default_well'.  If this is
%              false (default), prescribed well transmissibilities will be
%              kept (and updated based on the selected refinement).
%              Otherwise, transmissibilities will be replaced with default
%              values (i.e. to be computed from the grid).
%
% RETURNS:
%   deck - defined deck
%
% EXAMPLE:
%
% SEE ALSO:
%   refineGRDECL

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

   opt.default_well = false;
   opt = merge_options(opt, varargin{:});
   
   % Directly copy specific fields
   fields = {'SOLUTION','REGIONS','PROPS','RUNSPEC','SUMMARY', ...
             'SCHEDULE', 'UnhandledKeywords'};
   
   for f = fields
      f = f{:}; %#ok
      if isfield(deck_in, f)
         deck.(f)=deck_in.(f);
      end
   end

   % Fix RUNSPEC
   deck.RUNSPEC.DIMENS = deck_in.RUNSPEC.DIMENS.*dim;
   deck.RUNSPEC.cartDims = deck.RUNSPEC.DIMENS;
      
   % Expanding cell-based fields
   deck = expand_field(deck, dim, 'SOLUTION', {'EQLNUM', 'PRESSURE', 'SWAT', ...
                                               'SOIL', 'SGAS'});

   deck = expand_field(deck, dim, 'REGIONS', {'SATNUM','FIPNUM','PVTNUM', ...
                                              'ROCKNUM', 'IMBNUM', 'EQLNUM'});

   deck = expand_field(deck, dim, 'PROPS', {'SWL', 'SWCR', 'SGU', 'SGL', ...
                                            'SGCR', 'SOWCR', 'SOGCR', 'SWU', ...
                                            'ISWCR', 'ISGU', 'ISWL', 'ISWU', ...
                                            'ISGL', 'ISOGCR', 'ISOWCR', ...
                                            'ISGCR', 'SWATINIT'});

   % Refine the grid
   deck.GRID = refineGrdecl(deck_in.GRID, dim);

   % Refine schedule (i.e. wells)
   if isfield(deck, 'SCHEDULE') 
      deck.SCHEDULE = refine_wells(deck.SCHEDULE, dim, opt.default_well);
   end
end

% -------------------------------------------------------------------------

function deck = expand_field(deck, dim, gname, fname)
   cartDims_new = deck.RUNSPEC.cartDims;
   cartDims_orig = cartDims_new ./ dim;
   
   if isfield(deck, gname)
      for f = fname
         f = f{:}; %#ok
         if isfield(deck.(gname), f)

            A = reshape(deck.(gname).(f), cartDims_orig);
            tmp=A(ceil((1 : cartDims_new(1)) / dim(1)), ...
                  ceil((1 : cartDims_new(2)) / dim(2)), ...
                  ceil((1 : cartDims_new(3)) / dim(3)));
            deck.(gname).(f) = tmp(:);
         end
      end
   end
end

% -------------------------------------------------------------------------

function schedule = refine_wells(schedule, dim, default_well)
   for j = 1:numel(schedule.control)
      cur_ctrl = schedule.control(j);
      if ~isfield(cur_ctrl, 'WELSPECS')
         continue;
      end

      % Updating lateral (i,j) cell indices of well positions
      for i=1:2
         E = cur_ctrl.WELSPECS(:,2+i);
         schedule.control(j).WELSPECS(:,2+i) = ...
             cellfun(@(x) dim(i).*(x-1) + ceil(dim(i)/2), E, 'unif', false);
      end

      % Updating lateral cell indices (i, j) for individual well completions
      for i=1:2
         E = cur_ctrl.COMPDAT(:,1+i);
         schedule.control(j).COMPDAT(:,1+i) = ...
             cellfun(@(x) dim(i).*(x-1) + ceil(dim(i)/2), E, 'unif', false);
      end

      % Updating vertical cell indices (k) for individual well completions
      vw_ind = [cur_ctrl.COMPDAT{:,13}] == 'Z';  % indicator for vertical wells

      E = cur_ctrl.COMPDAT(:, 4); % start index 
      schedule.control(j).COMPDAT(vw_ind, 4) = ...
          cellfun(@(x) dim(3) .* (x-1) + 1, E(vw_ind), 'unif', false);
      schedule.control(j).COMPDAT(~vw_ind, 4) = ...
          cellfun(@(x) dim(3) .* (x-1) + ceil(dim(3)/2), E(~vw_ind), 'unif', false);

      E = cur_ctrl.COMPDAT(:, 5); % end index
      schedule.control(j).COMPDAT(vw_ind ,5) = ...
          cellfun(@(x) dim(3) .* x, E(vw_ind), 'unif', false);
      schedule.control(j).COMPDAT(~vw_ind, 5) = ...
          cellfun(@(x) dim(3) .* (x-1) + ceil(dim(3)/2), E(~vw_ind), 'unif', false);

      % Repeat perforations in the right direction and change fields
      compdat_old = schedule.control(j).COMPDAT;
      schedule.control(j).COMPDAT = [];

      % loop over perforations
      for i = 1:size(compdat_old, 1) 
         gg = compdat_old(i,:); % the current perforation
         tmp_str = struct('X', 1, 'Y', 2, 'Z', 3);
         mydim = tmp_str.(cell2mat(gg(13)));

         if mydim ~= 3 % z-direction already taken care of above

            % repeat the current completion once per refinement
            gg = repmat(gg, dim(mydim), 1); 
            gg(:, 1+mydim) = arrayfun(@(x) x, ...
                                      cell2mat(gg(:, 1+mydim)) + ...
                                      (0:dim(mydim)-1)', 'unif', false);
         end

         if (default_well)
            gg(:, [8, 10]) = {'-1'};
         else
            nz = [gg{:, 8}] > 0;
            if any(nz) % trans
               gg(nz, 8) = arrayfun(@(x) x, [gg{nz, 8}] / dim(mydim), 'unif', false); 
            end

            nz = [gg{:, 10}] > 0;
            if any(nz) % KH
               gg(nz, 10) = arrayfun(@(x) x, [gg{nz, 10}] / dim(mydim), 'unif', false); 
            end
         end

         schedule.control(j).COMPDAT=[schedule.control(j).COMPDAT; gg];

      end
   end
end
