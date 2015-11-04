%% 1D Antiform Caprock with Upscaling of Subscale Structures
% We consider a simple 1D antiform aquifer with a top surface given by the
% following expression
%      z = D - L1 * sin(x/L1)*tan(phi) + A sin(2*pi*x/L2)
% We study three different upscaling methods for the small-scale
% structures:
%  - An accretion layer
%  - An analytic approximation in terms of squares
%  - An analytic approximation in terms of a sinus wave
           
runStandardModel('data/upscalingExample1Data', @plotUpscalingFigs             , ...
                 'A', [2]                                                     , ...
                 'subscale_types', {'smooth', 'inf_rough', 'square', 'sinus'} , ...
                 'residual', true                                             , ...
                 'dis_types', {'none'});

             