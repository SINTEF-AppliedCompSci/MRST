classdef (InferiorClasses = {?DiagonalJacobian,?FixedWidthJacobian}) ConservationLawJacobian
    % Very experimental class for Jacobian of a conservation law
    properties
        flux
        accumulation
        divergenceOptions
    end
    
    methods
        function D = ConservationLawJacobian(acc, flux, div_options)
            D.accumulation = acc;
            D.flux = flux;
            D.divergenceOptions = div_options;
        end

        function s = sparse(D)
            s = discreteDivergenceDiagonalJac(D.accumulation, D.flux, D.divergenceOptions, false);
        end
        
        function [x, D] = diagMult(v, x, D)
            x = x.sparse();
            [x, D] = diagMult(v, x, D);
        end

        function D = mtimes(D, V)
            if isa(D, 'ConservationLawJacobian')
                D = D.sparse();
            else
                V = V.sparse();
            end
            D = mtimes(D, V);
        end
        
        function D = uminus(D)
            D.flux = -D.flux;
            D.accumulation = -D.accumulation;
        end

        function u = plus(u,v)
            if isa(u, 'ConservationLawJacobian')
                if isempty(u.accumulation)
                    u.accumulation = v;
                else
                    u.accumulation = u.accumulation + v;
                end
            else
                u = plus(v, u);
            end
        end
        
        
%         function varargout = matrixDims(D, n)
%             dims = [numel(D.mexPrelim.cellIndex), D.flux.dim(2)];
%             if nargout == 1
%                 varargout{1} = dims;
%                 if nargin > 1
%                     varargout{1} = varargout{1}(n);
%                 end
%             else
%                 varargout = {dims(1), dims(2)};
%             end
%         end

        function u = horzcat(varargin)
            u = ConservationLawJacobian.cat(2, varargin{:});
        end

        function u = vertcat(varargin)
            u = ConservationLawJacobian.cat(1, varargin{:});
        end
    end
    
    methods(Static)
        function out = cat(dim, varargin)
            varargin = cellfun(@sparse, varargin, 'UniformOutput', false);
            out = cat(dim, varargin{:});
        end
    end
end

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
