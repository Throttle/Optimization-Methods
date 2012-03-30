classdef optimizer < handle
    properties (Constant = true)
        delay = 0.3,
        point_size = 33
    end
    
    properties
        count,
        points
    end
    
    methods
        %% вычисляемая функция
        function y = f(this,x)
            this.count = this.count + 1;
            y = sinh((3*x^4 - x + sqrt(17) - 3) / 2) + sin((5^(1/3) * x^3 - 5^(1/3) * x + 1 - 2 * 5^(1/3))/(-x^3 + x + 2));
        end
        
        %% сброс счестчика вызова функции
        function reset(this)
            this.count = 0;
        end
        
        %% построение графика функции
        function draw(this, a, b, n)
            step = (b - a) / n;
            x = a : step : b;
            y = zeros(n + 1, 1);
            
            for i = 1 : n + 1
                y(i,1) = this.f(x(i));
            end
            figure(1);            
            plot(x, y(:,1));
        end
        
        %% прорисовка точки
        function draw_point(this, x, y)
            scatter(x, y, this.point_size, 'filled');
            pause(this.delay);
        end
                
        %% метод поразрядного поиска
        function [x, y] = RadixSearch(this, a, b, eps, k)
            fprintf('==== Метод поразрядного поиска ====\nОтрезок: a=%f, b=%f\nТочность: eps=%f;\nk=%d\n', [a,b,eps,k]);
            this.reset();
            delta = (b-a) / k;
            x = a;
            f0 = this.f(x);
            fprintf('Начало поиска: delta=%.7f; x=%.7f; f(x)=%.7f\n', [delta,x,f0]);
            while abs(delta) > eps
                while (x >= a) && (x <= b)
                    this.draw_point(x, f0);
                    x = x + delta;
                    f1 = this.f(x);
                    fprintf('f(%.7f)=%.7f\n', [x, f1]);
                    
                    if f1 > f0
                        break;
                    end

                    f0 = f1;
                end
                f0 = f1;
                delta = - delta / k;
                fprintf('Меняем направление поиска: delta=%.7f\n', delta);
            end
            disp('Поиск закончен: достигнута требуемая точность');
            y = f1;
            fprintf('xmin=%.7f\nf(xmin)=%.7f\n', [x,y])
            fprintf('Кол-во вызовов: %d\n', [this.count])
        end
        
        %% метод золотого сечения
        function [x, y, res_a, res_b] = GoldenSectionSearch(this, a, b, eps, stop, need_reset)
            if nargin < 5
                stop = -1;
            elseif nargin < 6
                this.reset();
            end
            
            fprintf('==== Метод золотого сечения ====\nОтрезок: a=%f, b=%f\nТочность: eps=%f;\n\tКол-во итераций: stop=%d.\n\n', [a,b,eps,stop]);
            %this.reset();
            tau = (sqrt(5) - 1) / 2;
            interval = b - a;
            
            x1 = a + (1 - tau) * interval;
            x2 = a + tau * interval;
            f1 = this.f(x1);
            f2 = this.f(x2);
            exec_count = 0;
            %fprintf('Начало поиска:\n\tотрезок: a=%.7f; b=%.7f;\n\tточки разбиения: f(%.7f)=%.7f; f(%.7f)=%.7f\n', [a, b, x1, f1, x2, f2]);
            while (abs(interval) > eps) && (exec_count - stop ~= 0)
                
                this.draw_point(x1, f1);
                this.draw_point(x2, f2);
                fprintf('-> Итерация %d\nИсходный отрезок:\n\ta=%.7f; b=%.7f;\n\tточки разбиения: f(%.7f)=%.7f; f(%.7f)=%.7f\n\n', [exec_count+1, a, b, x1, f1, x2, f2]);
                
                if f1 < f2
                    b = x2;
                    interval = b-a;
                    x2 = x1;
                    x1 = a + (1 - tau) * interval;
                    f2 = f1;
                    f1 = this.f(x1);
                else
                    a = x1;
                    interval = b - a;
                    x1 = x2;
                    x2 = a + tau * interval;
                    f1 = f2;
                    f2 = this.f(x2);
                end
                fprintf('Новый отрезок:\n\ta=%.7f; b=%.7f;\n\tзначения функции: f(a)=%.7f; f(b)=%.7f\n\n', [a, b, f1, f2]);
                exec_count = exec_count + 1;
            end
            disp('Поиск закончен: достигнута требуемая точность');
                        
            if f1 < f2
                x = x1;
                y = f1;
            else
                x = x2;
                y = f2;
            end
            
            fprintf('xmin=%.7f\nf(xmin)=%.7f\n', [x,y])
            fprintf('Кол-во вызовов: %d\n', [this.count])
            res_a = a;
            res_b = b;
            fprintf('=============\n\n')
        end
        
        %% Модифицированный метод Ньютона
        function [x, y] = Newton(this, a, b, eps, h)
            fprintf('==== Модифицированный метод Ньютона ====\nОтрезок: a=%f, b=%f\nТочность: eps=%f;\n', [a, b, eps]);
            this.reset();
            x = a;
            while 1
                f1 = this.f(x - h);
                f2 = this.f(x + h);
                f  = this.f(x);
                x_prev = x;
                fprintf('Точка: x=%.7f; f(x)=%.7f\n', [x, f]);
                this.draw_point(x, f);
                p = h * (f2 - f1) / (2 * (f2 - 2*f + f1));
                
                % модифицированный метод Ньютона
                alpha = 1;
                x_new = x - alpha * p;
                while (x_new < a) || (x_new > b)
                    fprintf('Точка (x=%.7f) за пределами отрезка [%.7f, %.7f]\n', [x_new, a, b]);
                    alpha = alpha / 2;
                    x_new = x - alpha * p;
                end
                x = x_new;
                
                if abs(x - x_prev) <= eps
                    break
                end
            end
            y = f;
            fprintf('xmin=%.7f\nf(xmin)=%.7f\n', [x,y])
            fprintf('Кол-во вызовов: %d\n', [this.count])
            fprintf('==== конец ====\n')
        end
        
        %%
        function result = r(this,x1,x2)
            result = x1^2 - x2^2;
        end
        
        %%
        function result = s(this,x1,x2)
            result = x1 - x2;
        end
        
        %% Метод парабол
        function [x, y] = Parabola(this, a, b, eps, golden_stop)
            if nargin < 5
                % отключаем выполнение метода золотого сечения
                golden_stop = 0;
            end
            fprintf('==== Метод парабол ====\nОтрезок: a=%f, b=%f\nТочность: eps=%f;\n', [a, b, eps]);
            this.reset();
            [xmin, fmin, a, b] = this.GoldenSectionSearch(a, b, eps, golden_stop, 0);
            
            x = [a, (a+b)/2, b];
            f = [this.f(x(1)), this.f(x(2)), this.f(x(3))];
            
            while abs(x(1) - x(3)) > eps
                x_dot = 0.5 * (...
                    f(1) * this.r(x(2), x(3)) + ...
                    f(2) * this.r(x(3), x(1)) + ...
                    f(3) * this.r(x(1), x(2)) ...
                    ) / (...
                    f(1) * this.s(x(2), x(3)) + ...
                    f(2) * this.s(x(3), x(1)) + ...
                    f(3) * this.s(x(1), x(2))...
                    );
                
                f_dot = this.f(x_dot);
                this.draw_point(x_dot, f_dot);
                
                if (abs(x(2) - x_dot) < eps)
                    break;
                end
                
                
                if (x_dot >= x(2)) && (x_dot <= x(3))
                    fprintf('Точка х*(%.7f) найдена внутри отрезка [x2,x3]([%.7f, %.7f]).\n', [x_dot, x(2), x(3)]);
                    if f_dot <= f(2)
                        fprintf('f(x*)=%.7f <= f(x2)=%.7f.\n', [f_dot, f(2)]);
                        x(1) = x(2);    f(1) = this.f(x(1));
                        x(2) = x_dot;   f(2) = this.f(x(2));
                    else
                        fprintf('f(x*)=%.7f > f(x2)=%.7f.\n', [f_dot, f(2)]);
                        x(3) = x_dot;   f(3) = this.f(x(3));
                    end
                elseif (x_dot >= x(1)) && (x_dot <= x(2))
                    fprintf('Точка х*(%.7f) найдена внутри отрезка [x1,x2]([%.7f, %.7f]).\n', [x_dot, x(1), x(2)]);
                    if f_dot <= f(2)
                        fprintf('f(x*)=%.7f <= f(x2)=%.7f.\n', [f_dot, f(2)]);
                        x(3) = x(2);    f(3) = this.f(x(3));
                        x(2) = x_dot;   f(2) = this.f(x(2));
                    else
                        fprintf('f(x*)=%.7f > f(x2)=%.7f.\n', [f_dot, f(2)]);
                        x(1) = x_dot;   f(1) = this.f(x(1));
                    end
                else
                    fprintf('Точка х*(%.7f) найдена вне отрезка [%.7f, %.7f]. Применяем метод золотого сечения.\n', [x_dot, x(1), x(3)]);
                    [xmin, fmin, a, b] = this.GoldenSectionSearch(x(1), x(3), eps, golden_stop, 0);
                    x = [a, (a+b)/2, b];
                    f = [this.f(x(1)), this.f(x(2)), this.f(x(3))];
                end
            end
            
            y = f_dot;
            x = x_dot;
            fprintf('xmin=%.7f\nf(xmin)=%.7f\n', [x,y])
            fprintf('Кол-во вызовов: %d\n', [this.count])
        end
    end
end