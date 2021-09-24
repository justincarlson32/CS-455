clc,clear
close all

qi = [1 5];

qj = [7,6];


sigmaNormResults = sigmaNormFun(qi-qj);

disp('For points [1 5] [7 6]');

disp('Sigma-Norm: ')
disp(sigmaNormFun(qi-qj));

disp('Eucl: ')
disp(norm(qi-qj));

disp('For points [3 8] [4 5]');

qi = [3 8];
qj = [4 5];

disp('Sigma-Norm: ')
disp(sigmaNormFun(qi-qj));

disp('Eucl: ')
disp(norm(qi-qj));

disp('For points [6 7] [5 5]');

qi = [6 7];
qj = [5 5];

disp('Sigma-Norm: ')
disp(sigmaNormFun(qi-qj));

disp('Eucl: ')
disp(norm(qi-qj));

function sigmaNormFun = sigmaNormFun(arg1)
    epsilon = 0.1; %(0,1)
    sigmaNormFun = (1/epsilon)*(sqrt(1 + epsilon*(norm(arg1,2)^2)) - 1);
end
