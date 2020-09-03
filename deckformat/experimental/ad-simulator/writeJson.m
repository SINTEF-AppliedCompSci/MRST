function writeJson(sopt,optfile)
fid = fopen(optfile,'wt');
fprintf(fid, '%s', jsonencode(sopt));
fclose(fid); 
end

