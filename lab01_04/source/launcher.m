clear;
clc;

optimizer = optimizer();

a = 0.0;
b = 1.0;

eps1 = 0.01;
eps2 = 0.0001;
eps3 = 0.000001;

optimizer.draw(a, b, 1000);
%[result, value] = optimizer.RadixSearch(a, b, eps3, 4);
%[result, value, res_a, res_b] = optimizer.GoldenSectionSearch(a, b, eps3, 0);
%[result, value] = optimizer.Newton(a, b, eps3, 1e-5);
[result, value] = optimizer.Parabola(a, b, eps3, 5);