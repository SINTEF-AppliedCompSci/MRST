%% 1D Antiform Caprock with Subscale Structures
% We consider a simple 1D antiform aquifer with a top surface given by the
% following expression
%      z = D - L1 * sin(x/L1)*tan(phi) + A sin(2*pi*x/L2)
% We study the 50 years of injection followed by 2000 years of migration
% for using two different geological models
%   - smooth caprock without small-scale traps (A=0)
%   - caprock with small-scale traps (A=2)
% and two different fluid models:
%   - plume migration without residual trapping
%   - plume migration with residual trapping

runStandardModel('data/residualExample1Data', @plotResidualFigs , ...
                 'A'             , [0 2]                        , ...
                 'residual'      , [false true]                 , ...
                 'dis_types'     , {'none'});