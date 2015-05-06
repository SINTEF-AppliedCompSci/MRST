function yesno = ask_user(question)

   yesno = [];
   
   while isempty(intersect(yesno, {'y', 'Y', 'n', 'N'}))
      yesno = input(question, 's');
   end
   yesno = strcmpi(yesno, 'y');
end
