classdef optimizer < handle
    properties (Constant = true)
        n = 2,
        delay = 0.1,
        point_size = 22,
        eps = 0.01,
        gl_H = 1e-4
    end
    
    properties
        count,
        points,
        color,
        m_rK
    end
    
    methods
        
        function y = h1(this, x)
            y = -2*x(1) +3* x(2) - 14;
        end

        function y = h2(this, x)
            y = -x(1) - 3;
        end

        function y = h3(this, x)
            y = (x(1) +2)^2 + (x(2) - 4)^2 - 15;
        end
        
        function y = f_b(this, x)
            y = this.f(x) - this.m_rK * this.bh(x);
        end
        
        function y = f_p(this, x)
            y = this.f(x) + this.m_rK * this.ph(x);
        end
        
        function y = bh(this, x)
           y = 1/this.h1(x) + 1/this.h2(x) + 1/this.h3(x);
        end
        
        function y = ph(this, x)
           f1 = this.h1(x);
           f2 = this.h2(x);
           f3 = this.h3(x);
           y = this.gPenalty(f1) + this.gPenalty(f2) + this.gPenalty(f3);
        end
        
        function g = gPenalty(this, fx)
            g = ((fx + abs(fx)) / 2) ^ 2;
        end
        
        %% вычисляемая функция
        function y = f(this, x)
            this.count = this.count + 1;
            % 1ая функция
            y = 4*x(1)*x(2) + 7*x(1)*x(1) + 4*x(2)*x(2) + 6*sqrt(5)*x(1) - 12*sqrt(5)*x(2) + 51;
        end
        
        %% сброс счестчика вызова функции
        function reset(this)
            this.count = 0;
        end
                
        %% построение графика функции в виде семейства линий уровня минимизируемых функций
        function draw(this, a1, b1, a2, b2, n, function_numb)
            step1 = (b1 - a1) / n;
            step2 = (b2 - a2) / n;
            x1 = a1 : step1 : b1;
            x2 = a2 : step2 : b2;
            %x2 = -10 : step : 2;
            f = zeros(length(x1), length(x2));
            
            for i = 1 : length(x1)
                for j = 1 : length(x2)
                    f(i,j) = this.f([x1(i), x2(j)]);
                end
            end
            
            for i = 1 : length(x1)
                for j = 1 : length(x2)
                    if function_numb == 1
                        bh(i,j) = this.ph([x1(i), x2(j)]);
                    else
                        bh(i,j) = this.bh([x1(i), x2(j)]);
                    end
                end
            end

            
            figure(1);
            levels = -30:5:50;
            [C,h] = contour(x2, x1, f, levels);
            set(h,'ShowText','on')
            levels = -10:2:10;
            [C,h] = contour(x2, x1, bh, levels);
            %[C,h] = contour(x2, x1, f);
            set(h,'ShowText','on')
            %colormap cool
        end
        
        %% прорисовка точки
        function draw_point(this, x, y)
            scatter(x, y, this.point_size, 'filled');
            pause(this.delay);
        end
        
        %% прорисовка линии
        function draw_line(this, p1, p2)
            line([p1(2), p2(2)], [p1(1), p2(1)]);
            pause(this.delay);
        end
        
        %% прорисовка симплекса
        function draw_simplex(this, S, color)
            if nargin < 3
                color = this.color;
            end
            line([S(2,2), S(1,2)], [S(2,1), S(1,1)], 'Color', color);
            line([S(3,2), S(1,2)], [S(3,1), S(1,1)], 'Color', color);
            line([S(3,2), S(2,2)], [S(3,1), S(2,1)], 'Color', color);
            pause(this.delay);
        end
        
       
        %% построение симлекса
        function S = get_simplex(this, x0, a)
            x = zeros(this.n+1, this.n);
            for i = 1 : this.n+1
                for j = 1 : this.n
                    if i == 1
                        x(i, j) = x0(j);
                    elseif i == (j + 1)
                        x(i, j) = x0(j) + a * (sqrt(this.n + 1) - 1)/(this.n * sqrt(2)); 
                    else
                        x(i, j) = x0(j) + a * (sqrt(this.n + 1) + this.n - 1)/(this.n * sqrt(2));
                    end
                end
            end
            S = x;
        end
                
        %% метод минимизации по правильному симплексу
        function [resx, resf] = simple_simplex(this, x0, a, eps)
            fprintf('==== Метод минимизации по правильному симплексу ====\nБазовая точка: x0=[%f, %f]\nТочность: eps=%f;\nДлина ребра: a=%d.\n\n', [x0(1),x0(2),eps,a]);
            this.color = 'b';
            this.reset();
            
            s = this.get_simplex(x0, a);
            this.draw_simplex(s);
                
            f = zeros(this.n+1, 1);    
            for i = 1 : this.n+1
                f(i) = this.f(s(i, :));
            end
            
            while true
                success = false;
                [sorted_f, unsorted_indexes] = sort(f);
                for k = length(unsorted_indexes) : -1 : 1
                    % вычисляем значение точки центра тяжести
                    maxindex = unsorted_indexes(k);
                    xc = sum(s) - s(maxindex,:);
                    this.draw_point(xc(2)/(this.n), xc(1)/(this.n));
                    new_point = [xc(2)*2/(this.n) - s(maxindex,2), xc(1)*2/(this.n) - s(maxindex,1)];
                    this.draw_point(new_point(1), new_point(2));
                    new_f = this.f([new_point(2), new_point(1)]);
                    if new_f > f(maxindex)
                        bad_simplex = s;
                        bad_simplex(maxindex,:) = [new_point(2), new_point(1)];
                        this.draw_simplex(bad_simplex, 'r');
                        continue;
                    else
                        f(maxindex) = new_f;
                        s(maxindex,:) = [new_point(2), new_point(1)];
                        this.draw_simplex(s);
                        success = true;
                        break;
                    end
                end
                
                if success == false
                    minindex = unsorted_indexes(1);
                    minpoint = s(minindex, :);
                    for i = 1:this.n+1
                        if i == minpoint
                            continue;
                        end
                        s(i,1) = minpoint(1) + this.delta * (s(i,1) - minpoint(1));
                        s(i,2) = minpoint(2) + this.delta * (s(i,2) - minpoint(2));
                    end
                    this.color = 'g';
                    this.draw_simplex(s);
                end
                
                summ = 0.0;
                m = sorted_f(1);
                for i = 1 : this.n
                    summ = summ + (f(i) - m)^2;
                end
                stopxx = sqrt((s(1,1) - s(2,1))^2 + (s(1,2) - s(2,2))^2) <= eps;
                stopxf = sqrt(summ / this.n) <= eps;
                
                if stopxf || stopxx
                    resx = s(unsorted_indexes(1),:);
                    resf = sorted_f(1);
                    fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n', [resx(1),resx(2),resf,this.count]);
                    break;
                end
            end
        end
        
        function w = grad(this, x)
            fx = this.f(x);
            w = [(fx - this.f([x(1) - this.gl_H, x(2)])) / this.gl_H;  (fx - this.f([x(1), x(2) - this.gl_H])) / this.gl_H];
        end
        
        function w = grad_h(this, x)
            fx = this.bh(x);
            w = [(fx - this.bh([x(1) - this.gl_H, x(2)])) / this.gl_H;  (fx - this.bh([x(1), x(2) - this.gl_H])) / this.gl_H];
        end
        
        function [resx, resf] = penalty(this, x0, eps)    
            this.reset()
            xK = x0;
            this.m_rK = 2;
            
            flag = true;
            this.draw_point(x0(2), x0(1));
                
            while flag 
                t = fminsearch(@this.f_p, xK, optimset('TolX', eps));
                this.draw_point(t(2), t(1));
                
                flag = abs(this.f(t) - this.f(xK)) > eps;
                this.m_rK = this.m_rK + 2;
                xK = t;
            end
            
            rx = xK;
            ry = this.f(xK);
            fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n', [rx(1),rx(2),ry,this.count]);
        end
        
        
        function [resx, resf] = barrier(this, x0, eps)
            this.reset()
            xK = x0;
            g = this.grad(xK);
            g_h = this.grad_h(xK);
            this.m_rK = sum(g .* g_h) / sum(g_h .^2);
            this.draw_point(x0(2),x0(1));
            flag = true;
            
            while flag
                [resx, resf] = fminsearch(@this.f_b, xK, optimset('TolX', eps));
                this.draw_point(resx(2),resx(1));
                flag = abs(this.f(resx) - this.f(xK)) > eps;
                this.m_rK = 0.5 * this.m_rK;
                xK = resx;
            end
            
            resx = xK;
            resf = this.f(xK);
            fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n', [resx(1),resx(2),resf,this.count]);
        end
       
        %% использование возможностей Optimization Toolbox Matlab
        %function [resx, resf] = fminsearch(this, x0, eps)
        %    this.reset();
        %    fprintf('==== Optimization Toolbox Matlab ====\nБазовая точка: x0=[%f, %f]\nТочность: eps=%f;\n\n', [x0(1),x0(2),eps]);
        %    [resx, resf] = fminsearch(@this.f, x0, optimset('TolX', eps));
        %    fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n\n', [resx(1),resx(2),resf,this.count]);
        %end
    end
end