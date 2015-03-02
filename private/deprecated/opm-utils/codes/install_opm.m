% matlab mess up the LD_PATH we need fix the LD path before in the system
% call the ldfix variable need to be set to the overide matlabs use on its
% one libstdc++.so.6 libgcc_s.so.1
ldfix='export LD_LIBRARY_PATH=/lib/x86_64-linux-gnu/:/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH;'
command=[ldfix, 'python install_opm.py']
a=system(command)
if(a>0)
   error(['Command failed: ', command])
end
