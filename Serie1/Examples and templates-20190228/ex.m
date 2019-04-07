%%It is a good idea to have a different function for each part of the homework.

%Unfortunately, MatLab does not allow M-files with more than one function callable
%from outside. To circumvent this problem and not to have too many M-files,
%one for each section of the homework,
%the following trick can be used.

%N.B. Remember that the name of the script must coincide with
%the name of the first function which is the only callable from outside
%(ex. from the console or from another script placed in the same folder)

function ex(n)
%External function choosing which section to provide the solution to.
%Call it from the console specifying the section to be executed.
%Ex: >> ex1(1)  to solve the first section,
%    >> ex1(2)  to solve the second section...

switch n
    case 1
        ex1_1()           
    case 2
        ex1_2()
    case 3
        ex1_3()
end
end


function ex1_1()

%Function solving the 1st section of the homework
    disp('Ciao')
end

function ex1_2()
%Function solving the 2nd section of the homework
    disp('Hello')
end

function ex1_3()
%Function solving the 3rd section of the homework
    disp('Salut')
end

%...As many functions as the number of homework sections...

% Moreover, you can have additional functions which are used only by other
% functions contained in this M-file.
