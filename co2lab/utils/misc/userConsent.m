function c = userConsent(question)
%userConsent Ask user for his/her consent to a Y/N question
%  c = userConsent(question) requests a text response from the user in the
%  form of a Y/N question. The function returns logical 1 (true) if the
%  the first letter in the user's response is a 'y' (or 'Y'), and returns
%  logical 0 otherwise. Default answer is 'Y'. 
%
%  Example:
%  userConsent('Do you want to continue') will produce the following prompt
%  "Do you want to continue? Y/N [Y]: "

answ = input([question '? Y/N [Y]: '], 's');
if isempty(answ)
    answ = 'Y';
end
c = strncmpi(answ, 'Y', 1);
end

