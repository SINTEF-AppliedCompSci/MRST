function [wellsols,states,reports,extra] = runDeckOPM(deckfile,varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

%mrstModule add ad-blackoil ad-core mrst-gui ad-props deckformat
opt=struct('outputdir','output',...
    'simulator','flow',...
    'force_timestep',true,...
    'read_states',true,...
    'use_ebos_style',false,...
    'do_adjoint',false,...
    'no_output',false,...
    'verbose',true,...
    'np',1,...
    'threads',2,...
    'lineartol',1e-4,...
    'strongdefaults',true,...
    'storagecache',true,...
    'linsolver',[],...
    'linearsolver',[],...
    'only_read',false,...
    'end_step',[]);
opt=merge_options(opt,varargin{:});
% for extra see .XXX.DEBUG of an opm run
if(isempty(opt.linearsolver))
   amg = struct('type','amg',...
                    'maxlevel',10,...
                    'coarsenTarget',4000,...
                    'smoother','ILU0', ...
                    'post_smooth', 1,...
                    'pre_smooth', 1,...
                    'gamma', 1,...
                    'alpha', 0.3,...
                    'beta', 1e-5,...
                    'relaxation', 1,...
                    'skip_isolated',0,...
                    'verbosity',10);
                
   coarsesolver=struct('solver','loopsolver','preconditioner',amg,'maxiter',1,'tol',1e-1);
   finesmoother=struct('type','ILU0','relaxation',1);
   cpr = struct('type','cpr',...
                 'weight_type','trueimpes',...
                'finesmoother',finesmoother,...
                'coarsesolver',coarsesolver,...
                'pressure_var_index',1,...
                'pre_smooth',1,...
                'post_smooth',1,...
                'verbosity',1);                
   linearsolver = struct('preconditioner',cpr,'solver','bicgstab','tol',1e-4,'verbosity',1); 
else  
  preconditioner=struct('type','ILU0','relaxation',1);
  linearsolver = struct('precontitioner',preconditioner,'solver','bicgstab','tol',1e-5,'maxiter',200,'verbosity',3); 
end

linearsolverfile = 'tmp_linearsolver.json';
writeJson(linearsolver,linearsolverfile);


if(false)
    linearsolverconf = [' --linear-solver-configuration=file ',...
        '--linear-solver-configuration-json-file=',linearsolverfile,' '];
else
    linearsolverconf=[];
end
if(opt.np >0)
    command = ['mpirun -np ',num2str(opt.np),' '];
else
    command = [];
end

if(opt.use_ebos_style)
    command=[command, opt.simulator,' --ecl-deck-file-name=',deckfile];
    %[mydir, case_name, ext]=fileparts(deckfile);
    [mydir, ~, ~]=fileparts(deckfile);
    opt.output=fullfile(mydir);
else
    if(opt.force_timestep)
        command=[command, opt.simulator,' --enable-well-operability-check=false --full-time-step-initially=false --enable-adaptive-time-stepping=false --flow-newton-max-iterations=20 '];
    else
        command=[command, opt.simulator,' --enable-well-operability-check=false --full-time-step-initially=false --enable-tuning=false '];
        % --time-step-control=iterationcount 
    end
    command = [command,' --solve-welleq-initially=true '];
    if(opt.storagecache)
        command =  [command,' --enable-storage-cache=true '];
    else
        command =  [command,' --enable-storage-cache=false '];
    end
    if(opt.strongdefaults)
        command = [command,' --use-inner-iterations-wells=false '];
        command = [command,' --tolerance-cnv=0.001 '];
        command = [command,' --tolerance-cnv-relaxed=0.01 '];
        command = [command,' --tolerance-wells=0.00001 '];
        command = [command,' --project-saturations=true '];
        command = [command,' --max-strict-iter=8 '];
        command = [command,' --ecl-enable-drift-compensation=false '];
        command = [command,' --relaxed-max-pv-fraction=0.0 '];
        command = [command,' --ds-max=0.2 '];
    end
    command = [command,linearsolverconf];
    
    command =[command,' --linear-solver-reduction=',num2str(opt.lineartol),' '];
    command=[command,' --threads-per-process=',num2str(opt.threads),' '];
    if(opt.do_adjoint)
       % needed for adjoint
       command = [command,' --solve-welleq-initially=true --use-adjoint=true'];
       %command = [command,' --matrix-add-well-contributions=true --preconditioner-add-well-contributions=true '];
       % adjust tolerances
       command = [command,...
			 ' --flow-linear-solver-verbosity=10  --linear-solver-max-iter=50 ',...
			 '--linear-solver-reduction=1e-10 ',... 
			  '--time-step-verbosity=10 ',...
                          '--tolerance-cnv=1e-6 ',...                        
			  '--tolerance-cnv-relaxed=1e-5 ',...                                                                    
			  '--tolerance-mb=1e-5 ',...                         
			  '--tolerance-pressure-ms-wells=1000 ',...
		          '--tolerance-well-control=1e-7 ',...
			 '--tolerance-wells=1e-4 '];
        % linear sovler
        %command=[command,'--use-umfpack=true ']
        %command=[command,'--use-amgcl=true --use-amgcl-drs=true ']
        %command=[command,'--use-amg=true -use-cpr=true -use-gmres=true ']
    else
      if(opt.verbose)
        command = [command,' --flow-linear-solver-verbosity=0 '];
      end
    end
    command = [command,' --output-dir=',opt.outputdir,' ',deckfile];
end
command =[command, opt.linsolver];
if(opt.no_output)
    command = [command,' >& /dev/null'];
end
if(not(opt.only_read))
    delete([fullfile(opt.outputdir,'extra_out'),'/*.txt'])
    disp('MRST runing flow')
    delete([fullfile(opt.outputdir),'*.UNSMRY']);
    if(~isempty(opt.end_step))
        command = [command,' --end-step=',num2str(opt.end_step)];
    end
    disp(command)
    a = system(command);
    if(a~=0)
        warning('funning flow failed to finnish')
        %delete([fullfile(opt.outputdir),'*.UNSMRY']);
    end
end
%% make normal output
if(nargout>0)
    opmdir=fullfile(opt.outputdir);
    [mydir,case_name,ext]=fileparts(deckfile);
    ofile=fullfile(opmdir,upper(case_name)); 
    %opm_smry = readEclipseSummaryUnFmt(ofile);
    [wellsols_smry, time_smry]  = convertSummaryToWellSols(ofile);
    ofile_cap=fullfile(opmdir,upper(case_name));
    casenm=fullfile(opt.outputdir,upper(case_name));
    if(opt.read_states)
        init = readEclipseOutputFileUnFmt([casenm,'.INIT']);
        grid = readEclipseOutputFileUnFmt([casenm,'.EGRID']);    
        [G, rock, N, T] = eclOut2mrst(init, grid);
        extra=struct('G',G,'rock',rock,'N',N,'T',T,'init',init, ...
                     'grid',grid,'syscomand', command);
        [states_opm,rstrt_opm] = convertRestartToStates(ofile_cap,G,...
        'use_opm',false,'includeWellSols',false,'wellSolsFromRestart',true,...
        'includeFluxes',false,'consistentWellSols',false);
        files=dir([fullfile(opt.outputdir,'extra_out'),'/*','txt']);
        sT=cell(numel(files),1);
        for i = 1:numel(files)
             states_opm{i}.sT = load(fullfile(files(i).folder,files(i).name));
        end
        extra.sT=sT;
    else
        states_opm=[];
        extra=[];
    end
    wellsols=wellsols_smry;
    reports = readInfoStep([ofile,'.INFOSTEP']);
    reports.ReservoirTime = time_smry;
    %reports=struct('ReservoirTime',time_smry,'adjoint',adjoint);
    states=states_opm;
end
