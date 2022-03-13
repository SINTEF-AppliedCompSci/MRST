function sample = generateRockSample(N, varargin)
% Generate a sample from random rock properties of dimension N
%
% SYNOPSIS:
%   sample = generateRockSample(N);
%   sample = generateRockSample(N, 'seed', seed, 'toVector', true);
%
% REQUIRED PARAMETERS:
%   N - Grid dimensions of the model
%
% OPTIONAL PARAMETERS:
%   'seed' - Seed for the random generator
%   'toVector' - Organize perm and poro as one-dimensional vectors.
%                Default: false
%   'avg_perm', 'std_perm' - mean and standard deviation of the log of the
%                            permeability. Default: 100*milli*darcy, 0.5
%   'avg_poro', 'std_poro' - mean and standard deviation of the porosity.
%                            Default: 0.5, 0.1
%   'min_poro' - minimum value of the porosity. Default: 1e-3
%   'covarianceFun' - Covariance function for the Gaussian process
%   'correlationLength' - Correlation length if no covarianceFun is
%                         provided. Default: 0.3
%   'smooth' - Smooth covariance in case no covarianceFun is provided.
%              Default: false
%
% RETURNS:
%   sample - a struct with fields 'perm' and 'poro', holding the grid's
%            permeability and porosity, respectively.
% 
% EXAMPLE:
%   ensembleExampleTutorial
%
% SEE ALSO
%   GaussianProcessND

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


    opt = struct('avg_perm'         , 100*milli*darcy, ...
                 'std_perm'         , 0.5            , ...
                 'avg_poro'         , 0.5            , ...
                 'std_poro'         , 0.1            , ...
                 'min_poro'         , 1e-3           , ...
                 'seed'             , []             , ...
                 'covarianceFun'    , []             , ...
                 'correlationLength', 0.3            , ...
                 'smooth'           , false          , ...
                 'toVector'         , false          );
    opt = merge_options(opt, varargin{:});

    fun = setCovarianceFun(opt);
    
    if ~isempty(opt.seed)
        rng(opt.seed);
    end
    
    if isstruct(N) && isfield(N, 'cells')
        % We have a grid instead of dimensions
        d = mean(N.cells.volumes.^(1/N.griddim));
        xmax = max(N.cells.centroids) + d/2;
        xmin = min(N.cells.centroids) - d/2;
        N = ceil((xmax - xmin)./d);
    end
    if numel(N) - nnz(N == 1) == 1
        p = GaussianProcess1D(N, fun);
    else
        p = GaussianProcessND(N, fun);
    end
    
    p = p - mean(p(:));
    p = p./std(p(:));
    
    perm = 10.^(log10(opt.avg_perm) + p.*opt.std_perm);
    poro = opt.avg_poro + p.*opt.std_poro;
    poro = max(poro, opt.min_poro);
    
    if opt.toVector
        perm = perm(:); poro = poro(:);
    end
    
    sample = struct('perm', perm, 'poro', poro);
end

function fun = setCovarianceFun(opt)
    fun = opt.covarianceFun;
    if isempty(opt.covarianceFun)
        if ~opt.smooth
            fun = @(x) exp(-sqrt(sum(x.^2,2))/opt.correlationLength);
        else
            fun = @(x) exp(-sum(x.^2,2)/opt.correlationLength);
        end
    end
end