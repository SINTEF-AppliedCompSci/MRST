function v = tbldispatch2(u, indstruct)
%
%
% SYNOPSIS:
%   function v = tbldispatch2(u, indstruct)
%
% DESCRIPTION: dispatch vector u according to indstruct (from tbl2 to pivot space)
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
    v = u(indstruct{2}.inds);
end
