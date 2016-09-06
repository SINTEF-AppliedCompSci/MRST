function grdecl = refineGrdecl(grdecl, dim, varargin)
%
% SYNOPSIS:
%   function grdecl = refineGrdecl(grdecl, dim, varargin)
%
% DESCRIPTION: function to refine a grdecl grid (Eclipse format)
%
%   - Refines the grid
%   - Correct cell data keywords,
%     keyword = {'ACTNUM', 'PERMX', 'PERMY', 'PERMZ', 'PORO'...
%     , 'MULTZ', 'MULTX', 'MULTY', 'EQLNUM'};
%   - Correct the FAULTS keyword, multipliers are not changed
%     which is only correct for refinement in z direction
%   - Add the keywords which do not change,
%     keyword = {'OIL', 'METRIC', 'WATER', 'EQUIL', 'ROCK', 'FLUID', 'START'}
%   - for wells (SHEDULING.control)
%       - change cell number in of wells
%       - NB !! do not remove perforations
%       - KH and well index not changed
%   - NB! not handle flowbased upscaling 'MULTZ','MULTX','MULTY'
%
%
% COMMENT : Function is not fully tested. To be used with care.
%
% PARAMETERS:
%   grdecl   - grid in eclipse format
%   dim      - [nx, ny, nz] defines ho much one should refine in each direction
%   varargin - optional parameters
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
% RETURNS:
%   grdecl - refined grid in eclipse format
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

    opt = struct('default_well', false);
    opt = merge_options(opt, varargin{:});
    [xyz, zcorn] = grdeclXYZ(grdecl);
    nx = dim(1); ny = dim(2); nz = dim(3);
    % if(nz>1)

    %% refining in z direction
    dz = zcorn(:, :, 2:2:end) - zcorn(:, :, 1:2:end);
    zcorn_new = zeros(size(zcorn, 1), size(zcorn, 2), size(zcorn, 3) * nz);
    if(nz>1)
        for i = 0:nz - 1
            zcorn_new(:, :, 2 * i + 1:2 * nz:end) = zcorn(:, :, 1:2:end) + i * dz / nz;
            zcorn_new(:, :, 2 * i + 2:2 * nz:end) = zcorn(:, :, 1:2:end) + (i + 1) * dz / nz;
        end
    else
        zcorn_new = zcorn;
    end

    %% refining in x direction
    xyz_new = xyz;
    if(nx>1)
        zcorn = zcorn_new;
        di = xyz(:, 2:end, :) - xyz(:, 1:end - 1, :);
        xyz_new = zeros(6, (size(xyz, 2) - 1) * nx + 1, size(xyz, 3));
        for i = 1:nx
            xyz_new(:, i:nx:end - 1, :) = xyz(:, 1:end - 1, :) + (i - 1) * di / nx;
        end
        xyz_new(:, end, :) = xyz(:, end, :);
        zcorn_new = zeros(size(zcorn, 1) * nx, size(zcorn, 2), size(zcorn, 3));
        dz = zcorn(2:2:end, :, :) - zcorn(1:2:end, :, :);
        if(nx>1)
            for i = 0:nx - 1
                zcorn_new(2 * i + 1:2 * nx:end, :, :) = zcorn(1:2:end, :, :) + i * dz / nx;
                zcorn_new(2 * i + 2:2 * nx:end, :, :) = zcorn(1:2:end, :, :) + (i + 1) * dz / nx;
            end
        else
            zcorn_new = zcorn;
        end
    end

    %% refining in y direction
    zcorn = zcorn_new;
    xyz = xyz_new;
    dj = xyz(:, :, 2:end) - xyz(:, :, 1:end - 1);
    if(ny>1)
        xyz_new = zeros(6, size(xyz, 2), (size(xyz, 3) - 1) * ny + 1);
        for i = 1:ny
            xyz_new(:, :, i:ny:end - 1) = xyz(:, :, 1:end - 1) + (i - 1) * dj / ny;
        end
        xyz_new(:, :, end) = xyz(:, :, end);
        zcorn_new = zeros(size(zcorn, 1), size(zcorn, 2) * ny, size(zcorn, 3));
        dz = zcorn(:, 2:2:end, :) - zcorn(:, 1:2:end, :);
        if(ny>1)
            for i = 0:ny - 1
                zcorn_new(:, 2 * i + 1:2 * ny:end, :) = zcorn(:, 1:2:end, :) + i * dz / ny;
                zcorn_new(:, 2 * i + 2:2 * ny:end, :) = zcorn(:, 1:2:end, :) + (i + 1) * dz / ny;
            end
        else
            zcorn_new = zcorn;
        end
    end

    grdecl_new.cartDims = grdecl.cartDims .* [nx, ny, nz];
    grdecl_new.COORD = xyz_new(:);
    grdecl_new.ZCORN = zcorn_new(:);
    grdecl_new.ACTNUM = zeros(grdecl.cartDims .* [nx, ny, nz]);

    cartDims = grdecl.cartDims;

    %% cell based fields
    keyword = {'ACTNUM', 'PERMX', 'PERMY', 'PERMZ', 'PORO'...
               , 'MULTZ', 'MULTX', 'MULTY', 'EQLNUM'};

    if ~isfield(grdecl, 'ACTNUM'),
        % Interpret no ACTNUM as all cells active.  This is consistent with
        % (e.g.) 'processGRDECL'.
        grdecl.ACTNUM = ones([prod(cartDims), 1]);
    end

    for i = 1:numel(keyword)
        if(isfield(grdecl, keyword{i}))
            A = reshape(grdecl.(keyword{i}), cartDims);
            grdecl_new.(keyword{i}) = A(ceil([1:nx * cartDims(1)] / nx), ceil([1:ny * cartDims(2)] / ny), ceil([1:nz * cartDims(3)] / nz));
            grdecl_new.(keyword{i}) = grdecl_new.(keyword{i})(:);
        end
    end

    %% not changed fields
    keyword = {'OIL', 'METRIC', 'WATER', 'EQUIL', 'ROCK', 'FLUID', 'START'};
    for i = 1:numel(keyword)
        if(isfield(grdecl, keyword{i}))
            grdecl_new.(keyword{i}) = grdecl.(keyword{i});
        end
    end

    %% fault
    if(isfield(grdecl, 'FAULTS'))
        faults = fields(grdecl.FAULTS);
        grdecl_new.FAULTS = grdecl.FAULTS;
        for i = 1:numel(faults)
            cells = grdecl_new.FAULTS.(faults{i}).cells;
            grdecl_new.FAULTS.(faults{i}).cells = (cells - 1) .* repmat(dim(ceil([1:2 * 3] / 2)), size(cells, 1), 1) + 1;
        end
    end

    if(isfield(grdecl, 'SCHEDULE'))
        grdecl_new.SCHEDULE = grdecl.SCHEDULE;
        for j = 1:numel(grdecl.SCHEDULE.control)
            if(isfield(grdecl.SCHEDULE.control(j), 'WELSPECS'))
                % wellspec mod
                for i = 1:2
                    E = grdecl.SCHEDULE.control(j).WELSPECS(:, 2 + i);
                    grdecl_new.SCHEDULE.control(j).WELSPECS(:, 2 + i) = cellfun(@(x) dim(i) .* (x - 1) + 1, E, 'unif', false);
                end
                for i = 1:4
                    E = grdecl.SCHEDULE.control(j).COMPDAT(:, 1 + i);
                    if(i<4)
                        fac = dim(i);
                    else
                        fac = dim(3);
                    end
                    grdecl_new.SCHEDULE.control.COMPDAT(:, 1 + i) = cellfun(@(x) fac .* (x - 1) + 1, E, 'unif', false);
                end
                %% repeat perforation in the right direction and change fields
                compdat_old = grdecl_new.SCHEDULE.control(j).COMPDAT;
                grdecl_new.SCHEDULE.control(j).COMPDAT = [];
                for i = 1:size(compdat_old, 1)
                    gg = compdat_old(i, :);
                    tmp_str = struct('X', 1, 'Y', 2, 'Z', 3);
                    mydim = tmp_str.(cell2mat(gg(end - 1)));
                    % mydim = dim(tmp_str.(cell2mat(gg(end - 1))));
                    gg = repmat(gg, dim(mydim), 1);
                    for d = 1:size(gg, 1)
                        if(mydim<3)
                            gg(d, 1 + mydim) = {cell2mat(gg(d, 1 + mydim)) + d - 1};
                            %%
                            if(mydim == 1)
                                gg(d, 3) = {cell2mat(gg(d, 3)) + floor(dim(2) / 2)};
                            else
                                gg(d, 2) = {cell2mat(gg(d, 2)) + floor(dim(1) / 2)};
                            end
                            gg(d, 4) = {cell2mat(gg(d, 4)) + floor(dim(3) / 2)};
                            gg(d, 5) = {cell2mat(gg(d, 5)) + floor(dim(3) / 2)};
                        else
                            gg(d, 4) = {cell2mat(gg(d, 4)) + d - 1};
                            gg(d, 5) = {cell2mat(gg(d, 5)) + d - 1};
                            % place the perforation in midle of cell
                            gg(d, 2) = {cell2mat(gg(d, 2)) + floor(dim(1) / 2)};
                            gg(d, 3) = {cell2mat(gg(d, 3)) + floor(dim(2) / 2)};
                        end
                        if(~opt.default_well)
                            if(cell2mat(gg(d, 8))>0)
                                gg(d, 8) = {cell2mat(gg(d, 8)) / mydim}; % trans
                            end
                            if(cell2mat(gg(d, 10))>0)
                                gg(d, 10) = {cell2mat(gg(d, 10)) / mydim}; % KH
                            end
                        else
                            gg(d, 8) = {' - 1'};
                            gg(d,10) = {' - 1'};
                        end
                    end
                    grdecl_new.SCHEDULE.control(j).COMPDAT = [grdecl_new.SCHEDULE.control(j).COMPDAT;gg];
                end
            end
        end
    end

    grdecl =grdecl_new;
end
