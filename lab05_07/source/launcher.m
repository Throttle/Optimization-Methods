clear;
clc;

optimizer = optimizer();

x0 = [0.0, sqrt(5)];
x02 = [3.0, 3.0];
a = 0.5;
a2 = 0.5;
eps1 = 0.01;
eps2 = 0.0001;
eps3 = 0.000001;
hold on;

%% Первая часть
%optimizer.draw(-6, 2, 0, 9, 100);
% использование метода минизации по правильному симплексу
%optimizer.simple_simplex(x0, a, eps3);

% использование метода минизации по деформируемому симплексу
%optimizer.downhill_simplex(x0, a, eps3);

% использование метода минимизациии случайного поиска в возвратом
%optimizer.returnrandsearch(x0, a, eps3);

%optimizer.fminsearch(x0, eps3);

%% Вторая часть
optimizer.draw(-2, 5, -2, 5, 100);

% использование метода минизации по правильному симплексу
%optimizer.simple_simplex(x02, a2, eps2);

% использование метода минизации по деформируемому симплексу
%optimizer.downhill_simplex(x02, a2, eps2);

% использование метода минимизациии случайного поиска в возвратом
%optimizer.returnrandsearch(x02, a2, eps2);

optimizer.fminsearch(x02, eps2);

hold off;