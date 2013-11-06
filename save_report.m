function save_report(mydir, report, num)
if(~isdir(mydir))
    mkdir(mydir);
end
myfile=sprintf('%s/report_%i',mydir,num);
save(myfile,'report');