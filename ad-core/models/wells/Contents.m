% WELLS
%
% Files
%   computeWellContributionsNew - Setup well (residual) equations and compute corresponding source terms.
%   setupWellControlEquations   - Setup well controll (residual) equations 
%   updateConnectionDP          - Explicit update of hydrostatic pressure difference between bottom hole
%   updateSwitchedControls      - Update active well controls based on well solution structure
%   updateSwitchedWellControls  - Check for violated well limits and switch controls. 
%   WellModel                   - Well model for three phase flow with black oil-style fluids

%{
#COPYRIGHT#
%}
