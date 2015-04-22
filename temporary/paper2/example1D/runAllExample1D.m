%% generate data
% loading mrst and CO2 lab should be ok?
%%{
% run for all cases
exampleVESlopingAquifer1D_runall_dis; 
% run linear case with upscaled model
exampleVESlopingAquifer1D_runall_smoothed; 
% do 1D top surface upscaling
simple1D_VE_upscaling.m
%}
clear all; 
% run all scripts for generating figures
scripts = {'plot1DGrid.m',...% ok fig 4 left
           'density1Dexample.m',...% ok fig 4 right        
           'example1D_grid_uscaledfluid.m',...% ok fig 11 right
           'example1D_makefig_cap.m',...% ok
           'example1D_makefig_dis_all',...% ok
           'example1D_makefig_fig1'...% ok
           'example1D_makefig_fig2'...% ok
           'example1D_makefig_fig2_b'...% ok
           'example1D_makefig_fig3'...% ok
           'example1D_plot_profile_ill',...
           'example1D_plot_profile_ill_reconstruct',...
           'plotVEReconstruction.m',...
           'example1D_plot_profile.m',...
           'plotVEFig19_20.m'
          }

for i = 1:numel(scripts)   
   run(scripts{i}) 
end