function [Gf, rockf] = setFacesToFractures(G, rock, faces, varargin)
%Set a collection of grid faces to fractures for use in discrete fracture models (DFM)
%
% SYNOPSIS:
%   [Gf, rockf] = setFacesToFractures(G, rock, faces)
%   [Gf, rockf] = setFacesToFractures(G, rock, faces, 'pn1', 'vn1', ...)
%
% DESCRIPTION:
%   Sets a subset of the grid faces to fractures for use with discrete
%   fracture modelling (DFM), and updates the rock accordingly.
%
% REQUIRED PARAMETERS:
%   G     - Grid structure with geometry
%
%   rock  - Rock data structure
%
%   faces - Vector of faces indices for faces to be set to fractures.
%           Either indices into (1:G.faces.num), or a logical mask.
%
% OPTIONAL PARAMETERS:
%   aperture - Aperture of fractures - should be scalar, or equal in length
%              to the 'faces' vector
%
%   poro     - Fracture porosity. Should be scalar, or equal in length to 
%              the 'faces' vector. Default: 0.5.
%
%   perm     - Fracture permeability. Should be empty, scalar or equal in
%              length to the 'faces' vector. If empty, permeability will be
%              computed from aperture.
%
%   'other'  - Any other fields that may be present in the rock structure
%              (such as thermal conductivity) can be assigned a
%              corresponding fracture value in the returned rock by optinal
%              imput values in the same way as aperture, poro, and perm.
%
% RETURNS:
%   Gf    - Grid with a a subset of the faces set to hybrid cells
%
%   rockf - Rock structure updated with extries for hybrid cells

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

    %---------------------------------------------------------------------%
    opt = struct('aperture'  , 1e-3, ...
                 'poro'      , 0.5 , ...
                 'perm'      , []  );
    [opt, faces] = processInput(G, rock, faces, opt, varargin{:});
    %---------------------------------------------------------------------%
    
    %---------------------------------------------------------------------%
    require dfm; % We use the dfm module for this
    % Flag faces that should be converted to fracture cells in the call to
    % 'addhybrid' below
    facetags        = zeros(G.faces.num, 1);
    facetags(faces) = 1;
    G.faces.tags    = facetags;
    tag = [];
    if isfield(G.faces, 'tag'), tag = G.faces.tag; end
    % Create hybrid grid
    aperture = opt.aperture(faces);
    Gf = addhybrid(G, faces, opt.aperture);
    if isfield(Gf.cells, 'indexMap')
        Gf.cells = rmfield(Gf.cells, 'indexMap');
    end
    % Set aperture
    Gf.cells.aperture = zeros(Gf.cells.num,1);
    Gf.cells.aperture(Gf.cells.hybrid > 0) = aperture;
    if ~isempty(tag)
        hf = zeros(Gf.faces.num - numel(tag) - nnz(faces),1);
        Gf.faces.tag = [tag; tag(faces); hf]; 
    end
    % Augment the rock structure
    rockf = augmentRock(rock, Gf, opt);
    %---------------------------------------------------------------------%

end

%-------------------------------------------------------------------------%
function [opt, faces] = processInput(G, rock, faces, opt, varargin)
    
    assert(size(rock.perm, 2) == 1); % @@only handle isotropic perm for now

    if ~islogical(faces)
        tmp        = false(G.faces.num, 1);
        tmp(faces) = true;
        faces      = tmp;
    end
    assert(numel(faces) == G.faces.num);
    assert(nnz(faces) > 0);

    [opt0, extra] = merge_options(opt, varargin{:});
    if isempty(opt0.perm)
        assert(~isempty(opt0.aperture))
        opt0.perm = opt0.aperture.^2/12;
    end
    
    % Handle optional extra fields
    for i = 1:2:numel(extra)
        opt0.(extra{i}) =  extra{i+1};
    end
    opt = struct();
    for name = fieldnames(opt0)'
        v = opt0.(name{1});
        if isscalar(v) || numel(v) == nnz(faces)
            opt.(name{1})        = zeros(G.faces.num,1);
            opt.(name{1})(faces) = v;
        else
            assert(numel(v) == G.faces.num);
            opt.(name{1}) = v;
        end
    end
    
end

%-------------------------------------------------------------------------%
function rock = augmentRock(rock, G, opt)

    % Identify hybrid cells (i.e. fracture cells)
    hybrid = G.cells.hybrid > 0;
    assert(size(rock.perm, 2) == 1); % @@only handle isotropic perm for now

    names = reshape(setdiff(fieldnames(opt), 'aperture'), 1, []);
    for name = names
        if ~isfield(rock, name{1})
            warning('`%s` is not a rock property - skipping this.', name{1}) 
            continue
        end        
        v = zeros(G.cells.num, 1);
        v(~hybrid) = rock.(name{1});
        v(hybrid)  = opt.(name{1})(G.cells.tags(hybrid));
        rock.(name{1}) = v;
    end
    
end
%-------------------------------------------------------------------------%
