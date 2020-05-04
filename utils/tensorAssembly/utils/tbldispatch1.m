function v = tbldispatch1(u, indstruct)
%
%
% SYNOPSIS:
%   function v = tbldispatch1(u, indstruct)
%
% DESCRIPTION: Low level function, dispatch vector u according to indstruct (from tbl1 to pivot space)
%
% PARAMETERS:
%   u         - vector to be dispatched
%   indstruct - structure that describes the dispatching
%
% RETURNS:
%   v - dispatched vector
%
% EXAMPLE:
%
% SEE ALSO: `crossIndexArray`, `IndexArrays`
%
    v = u(indstruct{1}.inds);
end
