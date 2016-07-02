function mrstNargInCheck(low, high, nargIn)
%Check number of input arguments to function
%
% SYNOPSIS:
%   mrstNargInCheck(low, high, nargIn)
%
% DESCRIPTION:
%   Fails execution (calls ERROR) unless actual number of input arguments
%   to calling function is between lower and upper limits inclusive.  This
%   function should usually not be called within a loop as it is
%   implemented in terms of DBSTACK.
%
% PARAMETERS:
%   low    - Minumum number of input arguments needed by calling function.
%            If empty (i.e., if ISEMPTY(low)), this limit is not checked.
%
%   high   - Maximum number of input arguments allowed by calling function.
%            If empty (i.e., if ISEMPTY(high)), this limit is not checked.
%
%   nargIn - Actual number of input arguments.  Must be scalar, integral
%            and non-negative.  Should be NARGIN unless there are special
%            circumstances.
%
% SEE ALSO:
%   nargin, nargchk, narginchk, error, dbstack.

%{
#COPYRIGHT#
%}

   assert (isnumeric(nargIn) && (numel(nargIn) == 1) && ...
           (mod(nargIn, 1) == 0) && (nargIn >= 0), ...
           'Argument count must be scalar, non-negative integer');

   st = dbstack;
   [caller, caller] = fileparts(st(2).file);                    %#ok<ASGLU>

   if ~isempty(low) && (nargIn < low),
      error('ArgCount:Low', ...
           ['Too few input arguments to function ''%''\n', ...
            'Expected at least %d arguments, but got %d'], ...
            caller, low, nargIn);
   end

   if ~isempty(high) && (nargIn > high),
      error('ArgCount:High', ...
           ['Too many input arguments to function ''%''\n', ...
            'Expected at most %d arguments, but got %d'], ...
            caller, high, nargIn);
   end   
end
