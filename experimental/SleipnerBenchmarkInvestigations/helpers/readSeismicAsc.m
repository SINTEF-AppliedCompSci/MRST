function data = readSeismicAsc(file)
%function to read the asci files for seismics. Assums a lot but is probably
%correct for this limited numer of files and give out put in resonable
%format
a = importdata(file,' ',20);
% cartdim linje
%%
lin=a.textdata{14}
b=strsplit(lin,':')
assert(strcmp(b{1},'# Grid_size'));
tmp=strsplit(b{end},'x')
assert(numel(tmp)==2);
grid_size=nan(1,2);
for i=1:2
    grid_size(i)=str2num(tmp{i});
end
lin=a.textdata{15}
b=strsplit(lin,':')
assert(strcmp(b{1},'# Grid_space'))
tmp=strsplit(b{end},',')
assert(numel(tmp)==4);
grid_space=nan(4,1);
for i=1:numel(tmp)
    grid_space(i)=str2num(tmp{i});
end
val=cell(3,1);
for i=1:numel(val)
   val{i} =nan(grid_size);
   ind=sub2ind(size(val{i}),a.data(:,4),a.data(:,5));
   val{i}(ind)=a.data(:,i);
end
% make cartesian
 x=linspace(grid_space(1),grid_space(2),grid_size(1));
 y=linspace(grid_space(3),grid_space(4),grid_size(2));
 
[X,Y]=meshgrid(x,y);
ind=isfinite(val{3}');
Xtmp=val{1}';
assert(all(X(ind)==Xtmp(ind)));
Ytmp=val{2}';
assert(all(Y(ind)==Ytmp(ind)));
data=struct('val',{val},'grid_size',grid_size,'grid_space',grid_space,'X',X,'Y',Y,'Z',val{3}');
end

