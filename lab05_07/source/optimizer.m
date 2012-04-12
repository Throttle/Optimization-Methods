classdef optimizer < handle
    properties (Constant = true)
        n = 2,
        delay = 0.3,
        alpha = 1,
        betta = 0.5,
        gamma = 2,
        delta = 0.5,
        point_size = 22,
        max_it = 500, 
        div_val = 10
    end
    
    properties
        count,
        points,
        color
    end
    
    methods
        %% вычисляемая функция
        function y = f(this, x)
            this.count = this.count + 1;
            % 1ая функция
            y = 4*x(1)*x(2) + 7*x(1)*x(1) + 4*x(2)*x(2) + 6*sqrt(5)*x(1) - 12*sqrt(5)*x(2) + 51;
            % 2ая функция
            %y = x(2)*x(2)*x(2) + 2 * x(1) * x(2) + 1 / sqrt(x(1)*x(2));
        end
        
        %% сброс счестчика вызова функции
        function reset(this)
            this.count = 0;
        end
        
        %% построение графика функции в виде семейства линий уровня минимизируемых функций
        function draw(this, a1, b1, a2, b2, n)
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
            
            figure(1);
            levels = -30:5:40;
            [C,h] = contour(x2, x1, f, levels);
            %figure(2);
            %surface(x2,x1,f);
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
        
        %% метод минимизации по деформируемому симплексу
        function [resx, resf] = downhill_simplex(this, x0, a, eps)
            fprintf('==== Метод минимизации по деформируемому симплексу ====\nБазовая точка: x0=[%f, %f]\nТочность: eps=%f;\nДлина ребра: a=%d.\n\n', [x0(1),x0(2),eps,a]);
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
                maxindex = unsorted_indexes(length(unsorted_indexes));
                f_h = sorted_f(this.n+1);
                f_g = sorted_f(this.n);
                f_l = sorted_f(1);
                
                x_h = s(unsorted_indexes(this.n+1),:);
                
                % вычисление центра тяжести симплекса
                xc = sum(s) - s(maxindex,:);
                this.draw_point(xc(2)/this.n, xc(1)/this.n);
                
                % отражение
                x_r = [(1 + this.alpha) * xc(2)/this.n - this.alpha * s(maxindex, 2), (1 + this.alpha) * xc(1)/this.n - this.alpha * s(maxindex, 1)];
                this.draw_point(x_r(1), x_r(2));
                f_r = this.f([x_r(2), x_r(1)]);
                
                if f_r <= f_l
                    x_e = [(1 + this.gamma * this.alpha) * xc(2)/this.n - this.gamma * this.alpha * s(maxindex, 2), (1 + this.gamma * this.alpha) * xc(1)/this.n - this.gamma * this.alpha * s(maxindex, 1)];
                    this.draw_point(x_e(1), x_e(2));
                    f_e = this.f([x_e(2), x_e(1)]);
                    
                    if f_e < f_l
                        f(maxindex) = f_e;
                        s(maxindex,:) = [x_e(2), x_e(1)];
                        this.draw_simplex(s);
                    else
                        %bad_simplex = s;
                        %bad_simplex(maxindex,:) = [new_point(2), new_point(1)];
                        %this.draw_simplex(bad_simplex, 'r');
                        f(maxindex) = f_r;
                        s(maxindex,:) = [x_r(2), x_r(1)];
                        this.draw_simplex(s);
                    end
                elseif (f_l < f_r) && (f_r <= f_g)    
                    f(maxindex) = f_r;
                    s(maxindex,:) = [x_r(2), x_r(1)];
                    this.draw_simplex(s);
                else
                    flag = false;
                    %if f_r < sorted_f(this.n+1)
                    if (f_h >= f_r) && (f_r > f_g)
                        flag = true;
                        s(maxindex,:) = [x_r(2), x_r(1)];
                        f(maxindex) = f_r;
                        x_h = [x_r(2), x_r(1)];
                        f_h = f_r;
                    end
                    if (f_r > f_h)
                        flag = true;
                    end
                    
                    if flag
                        % сжатие
                        % вычисление центра тяжести симплекса
                        %xc = sum(s) - s(maxindex,:);
                        %this.draw_point(xc(2)/this.n, xc(1)/this.n);
                        x_s = [(1 - this.betta) * xc(2)/this.n + this.betta * x_h(2), (1 - this.betta) * xc(1)/this.n + this.betta * x_h(1)];
                        f_s = this.f([x_s(2), x_s(1)]);
                        this.draw_point(x_s(1), x_s(2));
                        if f_s < f_h
                            s(maxindex,:) = [x_s(2), x_s(1)];
                            f(maxindex) = f_s;
                            x_h = [x_s(2), x_s(1)];
                            f_h = f_s;
                        else
                            minindex = unsorted_indexes(1);
                            minpoint = s(minindex, :);
                            for i = 1:this.n+1
                                if i == minindex
                                    continue
                                end
                                s(i,1) = minpoint(1) + this.delta * (s(i,1) - minpoint(1));
                                s(i,2) = minpoint(2) + this.delta * (s(i,2) - minpoint(2));
                            end
                        end
                    end
                end
                this.draw_simplex(s);
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
        
        function [x, y] = returnrandsearch(this, x0, a, eps)
            fprintf('==== Метод минимизации случайного поиска с возвратом при неудачном шаге ====\nБазовая точка: x0=[%f, %f]\nТочность: eps=%f;\n\n', [x0(1),x0(2),eps]);
            this.count = 0;
            
            lastpoint = x0;
            
            xk = x0;
            l = a;
            
            k = 0;
            j = 1;
            fk = this.f(xk);
            
            this.draw_point(xk(2), xk(1));
            lastpoint=xk;
                        
            while true
                r = this.rand();
                xNew = xk + l * r / sqrt(sum(r.*r));
                fNew = this.f(xNew);
                
                if (fNew < fk)
                    xk = xNew;
                    fk = fNew;
                    k = k + 1;
                   
                    this.draw_point(xk(2), xk(1));
                    this.draw_line(lastpoint, xk);
                    lastpoint=xk;
            
            
                    if (k > this.max_it) %Максимальное число шагов
                        x = xk;
                        y = fk;
                        fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n', [x(1),x(2), y,this.count]);
                        return;
                    end
                    j = 0;
                else
                    if (j < this.div_val)
                        j = j + 1;
                    else
                        if (l < eps) %Превышено количество неудачных попыток с заданным l - уменьшаем l 
                            x = xk;
                            y = fk;
                            fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n', [x(1),x(2), y,this.count]);
                            return;
                        else
                            l = l * this.delta;
                            j = 1;
                        end
                    end
                end
            end
			
			
                    
        end
        
        function [r] = rand(this)
            r = 2*randn(1, this.n) - ones(1, this.n);
        end
       
        %% использование возможностей Optimization Toolbox Matlab
        function [resx, resf] = fminsearch(this, x0, eps)
            fprintf('==== Optimization Toolbox Matlab ====\nБазовая точка: x0=[%f, %f]\nТочность: eps=%f;\n\n', [x0(1),x0(2),eps]);
            [resx, resf] = fminsearch(@this.f, x0, optimset('TolX', eps));
            fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n', [resx(1),resx(2),resf,this.count]);
        end
    end
end