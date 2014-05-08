%% Basic capacity estimates for the CO2 atlas data
% In this test, we construct top-surface grids for all the parts of the
% atlas data set that we are able to do this for. Altogether, this should
% give seventeen different grids. We then compeare the result with a saved results.
% It is bases on basicCapasityAtlas from paper ??
try
   require coarsegrid deckformat mex
catch %#ok<CTCH>
   mrstModule add coarsegrid deckformat mex
end


%% Compute trapping capacity for all aquifers
% To better report progress, we first load a low resolution version to get
% names of all aquifers. Then we load and process the full-resolution
% versions using both the cell-based and node-based methods.
profile off
profile on 
grdecls = getAtlasGrid('coarsening',10);
ng = numel(grdecls);
res = cell(ng,1);

for i=1:ng
    %%
    %profile off
    %profile on  
   fprintf('------------------------------------------------\n');
   fprintf('Processing %s ....\n', grdecls{i}.name);
   grdecl  = getAtlasGrid(grdecls{i}.name, 'coarsening', 1);
   G       = mprocessGRDECL(grdecl{1});
   G       = mcomputeGeometry(G(1));
   Gt      = topSurfaceGrid(G);
   tan     = trapAnalysis(Gt, false);
   tac     = trapAnalysis(Gt, true);
   
   res{i}.name      = grdecls{i}.name;
   res{i}.cells     = Gt.cells.num;
   res{i}.zmin      = min(Gt.cells.z);
   res{i}.zmax      = max(Gt.cells.z);
   res{i}.volume    = sum(G.cells.volumes);
   res{i}.ctrapvols = volumesOfTraps(Gt,tac);
   res{i}.ccapacity = sum(res{i}.ctrapvols);
   res{i}.ntrapvols = volumesOfTraps(Gt,tan);
   res{i}.ncapacity = sum(res{i}.ntrapvols);
   fprintf('done\n');
   %profile off;profile viewer
end
%%
profile off;profile viewer
%%
filename=mfilename('fullpath');
mydir = fileparts(filename);
res_saved=load(fullfile(mydir,'data/capasity_estimate_hmn.mat'));
res_saved=res_saved.res;
myfields=fieldnames(res_saved{1});
ok=true;
for i =1:numel(myfields)
   data=res{i}.(myfields{i});            
   if( isnumeric(data))
     ok = (res_saved{i}.(myfields{i})-data)/data < 1e-10;
     %assert((res_saved{i}.(myfields{i})-data)/data < 1e-2)
   end
end
if(ok)
    disp('Capasity estimate is consistent with previous results')
end