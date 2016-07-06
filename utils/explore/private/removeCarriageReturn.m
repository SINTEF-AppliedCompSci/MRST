function s = removeCarriageReturn(s)
%Internal function (undocumented)

%{
#COPYRIGHT#
%}

   if ispc,
      s = regexprep(s, '\r\n', '\n');
   end
end
