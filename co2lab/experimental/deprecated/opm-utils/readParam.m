function s=readParam(myfile)
myfile=fileread(myfile);
k = regexp(myfile,'([./\w]+)=([-./\w]+)\n','tokens');
k = vertcat(k{:}) .';
for i=1:size(k,2)
   b=k{1,i};
   b(find(b=='/'))='_';   
   if(b(1)=='_')
      b=b(2:end);
   end
   k{1,i}=b;
end   
s = struct(k{:});
