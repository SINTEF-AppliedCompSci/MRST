function nsteps = numberOfsteps(output_dir)
   press_files = dir(fullfile(output_dir, 'pressure', '*.txt'));
   nsteps      = numel(press_files);
end

