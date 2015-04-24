function varargout = evalMultipleRegions(fun, reg, x, varargin)
%Evaluate fluid functions in multiple regions
%
% SYNOPSIS:
%   output = evalMultipleRegions(fun, region, x)
%   output = evalMultipleRegions(fun, region, x, ...)
%
% PARAMETERS:
%   fun - Cell array, one element for each region, of function handles.
%         Each function handle, fun{i}, is responsible for computing a
%         specific fluid quantity (e.g., viscosity or relative
%         permeability) in the corresponding region, 'i'.
%
%         Each evaluation function is assumed to support a calling
%         interface of the form
%
%             [y{1:n}] = fun{i}(x)
%             [y{1:n}] = fun{i}(x, ...)
%
%         Where 'n' is an unspecified number of output values.  This
%         function implicitly uses nargout bumping.  See
%         <url:http://blogs.mathworks.com/loren/2009/04/14/
%          convenient-nargout-behavior/> for details.
%
%   reg - Region numbers.  One positive, integral number in 1:NUMEL(fun),
%         for each cell in the discretised reservoir model.  In the case of
%         relative permeability function evaluation, this is the SATNUM or
%         IMBNUM data.
%
%   x   - Evaluation variable.  If numeric, (e.g., saturations), then
%         values restricted to each region will be passed into the
%         corresponding, individual evaluation functions.  In this case,
%         the 'x' values must be ordered such that a row corresponds to a
%         specific grid cell.
%
%         If 'x' is a structure (e.g., a 'state' structure as defined by
%         function 'initResSol'), then the variable will be passed
%         unmodified into each evaluation function.
%
%   ... - Other parameters.  OPTIONAL.  These will be passed untouched into
%         the evaluation functions.  The evaluation functions must change
%         their behaviour accordingly.
%
% RETURNS:
%   output - Return values.  Each of the 'n' return values ('n' determined
%            by caller) is a numeric array with a number of rows
%            corresponding to the total number of cells in the discretised
%            reservoir geometry.  The number of columns of each output is
%            determined by the corresponding output of the individual
%            evaluation routines.
%
% SEE ALSO:
%   fluid_structure.

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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


   nout = max(nargout, 1);
   nreg = max(reg);

   varargout = cell([1, nout]);
   t         = cell([1, nout]);

   for r = 1 : nreg,
      i  = reg == r;   % May be too expensive.  Profile!
      nc = sum(i);

      if isnumeric(x),
         [t{:}] = fun{r}(x(i,:), varargin{:});
      else
         [t{:}] = fun{r}(x     , varargin{:});
      end

      for out = 1 : nout,
         assert (isnumeric(t{out}), ...
                 'Evaluation routines must produce numeric output only.');

         if size(t{out}, 1) == 1,
            t{out} = repmat(t{out}, [nc, 1]);
         end

         if isempty(varargout{out}),
            varargout{out} = zeros([numel(reg), size(t{out}, 2)]);
         end

         assert (size(t{out}, 1) == nc, ...
                ['Evaluation routines must produce output restricted ', ...
                 'to region only.']);

         varargout{out}(i,:) = t{out};
      end
   end
end
