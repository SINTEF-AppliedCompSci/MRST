function value = getOpmField(output_dir,myfield,step)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
nsteps=numberOfsteps(output_dir);
assert(nsteps<step);
myfile=sprintf('%s-%03d.dat',[dir_name,filesep,myfield],nn);
value=load(myfile);
end

