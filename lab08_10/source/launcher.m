clear;
clc;

optimizer = optimizer();

x01 = [0.0, sqrt(5)];
x02 = [3.0, 3.0];
eps1 = 0.01;
eps2 = 0.0001;
eps3 = 0.000001;
hold on;

%% Первая часть
optimizer.draw(100);
optimizer.conjugate_gradient(x02, eps3);
%optimizer.newton_with_approximation(x01, eps1);
%optimizer.DFP(x01, eps1);
%optimizer.DFP(x01, eps2);
%optimizer.DFP(x01, eps3);
%optimizer.fminsearch(x02, eps1);

%% Вторая часть
%optimizer.draw(-2, 5, -2, 5, 100);

% использование метода минизации по правильному симплексу
%optimizer.simple_simplex(x02, a2, eps3);

% использование метода минизации по деформируемому симплексу
%optimizer.downhill_simplex(x02, a2, eps3);

% использование метода минимизациии случайного поиска в возвратом
%optimizer.rand_search(x02, a2, eps3);
%
%optimizer.fminsearch(x02, eps3);

hold off;