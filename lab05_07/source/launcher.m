clear;
%clc;

optimizer = optimizer();

x0 = [0.0, sqrt(5)];
a = 0.4;
eps1 = 0.01;
eps2 = 0.0001;
eps3 = 0.000001;
hold on;

%% ������ �����
optimizer.draw(-6, 2, 0, 9, 100);
% ������������� ������ ��������� �� ����������� ���������
%optimizer.simple_simplex(x0, a, eps3);

% ������������� ������ ��������� �� �������������� ���������
%optimizer.downhill_simplex(x0, a, eps3);

%
optimizer.returnrandsearch(x0, a, eps3);

%optimizer.fminsearch(x0, eps3);
hold off;