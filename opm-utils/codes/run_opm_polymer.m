function output_dir=run_opm_polymer(param,executable)
codes_polymer={'polymer_reorder','sim_poly2p_incomp_reorder'};
if(nargin==2)
   assert(any(strcmp(executable,codes_polymer)));
elseif(nargin==1)
   executable=codes_polymer{1};
else
   param=[];
   executable=codes_polymer{1};
end
if(~isempty(param))
   input_filename=fullfile('test_run.param');
   paramStructToParamFile(param,input_filename);
end
%%
% matlab mess up the LD_PATH we need fix the LD path before in the system
% call the ldfix variable need to be set to the overide matlabs use on its
% one libstdc++.so.6 libgcc_s.so.1
ldfix='export LD_LIBRARY_PATH=/lib/x86_64-linux-gnu/:/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH;'
if(~isempty(param))
   command=[ldfix, fullfile(opm_code_dir(),'builds/opt/opm-polymer/examples/',executable),' ',input_filename];
else
   command=[ldfix, fullfile(opm_code_dir(),'builds/opt/opm-polymer/examples/',executable)];
end
a=system(command)
if(a>0)
   error(['Command failed: ', command])
else
   if(isfield(param,'output_dir'))
      output_dir=param.output_dir;
   else
      output_dir='output';
   end
   %param_out=paramToStruct(fullfile(output_dir,'spu_2p.param'));
end
