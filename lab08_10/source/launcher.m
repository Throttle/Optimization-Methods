clear;
clc;

optimizer = optimizer();

x01 = [0.0, sqrt(5)];
x02 = [3.0, 3.0];
eps1 = 0.01;
eps2 = 0.0001;
eps3 = 0.000001;
hold on;

%% ������ �����
optimizer.draw(100);
optimizer.conjugate_gradient(x02, eps3);
%optimizer.newton_with_approximation(x01, eps1);
%optimizer.DFP(x01, eps1);
%optimizer.DFP(x01, eps2);
%optimizer.DFP(x01, eps3);
%optimizer.fminsearch(x02, eps1);

%% ������ �����
%optimizer.draw(-2, 5, -2, 5, 100);

% ������������� ������ ��������� �� ����������� ���������
%optimizer.simple_simplex(x02, a2, eps3);

% ������������� ������ ��������� �� �������������� ���������
%optimizer.downhill_simplex(x02, a2, eps3);

% ������������� ������ ������������ ���������� ������ � ���������
%optimizer.rand_search(x02, a2, eps3);
%
%optimizer.fminsearch(x02, eps3);

hold off;