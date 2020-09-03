function json = readJSon(fname)
fid = fopen(fname);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
%%
json = jsondecode(str);
end

