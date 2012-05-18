clear;
clc;

optimizer = optimizer();

x0 = [0.0, sqrt(5)];
x02 = [-1.0, 8.0];
a = 0.5;
a2 = 0.5;
eps1 = 0.001;
hold on;

%% Первая часть
function_numb=2;
optimizer.draw(-6, 2, 0, 9, 100,function_numb);
% использование метода минизации по правильному симплексу
if function_numb == 1
    optimizer.penalty(x02, eps1);
else
    %optimizer.barrier(x0, eps1);
    optimizer.barrier([-1.0, 1.0], eps1);
    %optimizer.barrier([0.0, 3.0], eps1);
end
%

% использование метода минизации по деформируемому симплексу
%optimizer.downhill_simplex(x0, a, eps3);

% использование метода минимизациии случайного поиска в возвратом
%optimizer.rand_search(x0, a, eps3);

%optimizer.fminsearch(x0, eps3);

hold off;