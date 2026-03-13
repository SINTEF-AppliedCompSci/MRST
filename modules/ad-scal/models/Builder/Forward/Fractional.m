function model = Fractional(model)
%
% DESCRIPTION: calculates the fractional flow of the model
%
% SYNOPSIS:
%   model = Fractional(model)
%
% PARAMETERS:
%   model - struct containing following fields:
%   - experiment
%   - satfun
%
% RETURNS:
%   model - struct containing following fields:
%   - satfun.fw: the calculated fractional flow
%
% ----------------------------------
% (c) 2020-2022
% Siroos Azizmohammadi
% Omidreza Amrollahinasab
% Montanuniversität Leoben, Austria
% Chair of Reservoir Engineering
% https://dpe.ac.at/
% ----------------------------------
%
%%
muW = model.experiment.fluid.muW;
muO = model.experiment.fluid.muNW;
s = model.satfun.s;
kr = model.satfun.kr;
for i = 1 : length(s)
    krs = kr(s(i));
    krw = krs(1); kro = krs(2);
    M  = (krw / kro) * (muO / muW);
    fw(i) = 1 / (1 + 1 / M);
end
model.satfun.fw = fw';