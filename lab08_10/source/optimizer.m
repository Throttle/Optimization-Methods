classdef optimizer < handle
    properties (Constant = true)
        n = 2,
        delay = 0.0,
        betta = 0.5,
        gamma = 2,
        delta = 0.5,
        point_size = 22,
        max_it = 500, 
        div_val = 10,
        function_numb = 2,
        psi_eps = 1e-10,
        eps = 0.1,
        conjugate_reset = 4,
    end
    
    properties
        count,
        points,
        color,
        phi_xk,
        phi_pk,
        alpha
    end
    
    methods
        
        
        %% вычисляемая функция
        function y = f(this, x)
            this.count = this.count + 1;
            if this.function_numb == 1
                % 1ая функция
                y = 4*x(1)*x(2) + 7*x(1)*x(1) + 4*x(2)*x(2) + 6*sqrt(5)*x(1) - 12*sqrt(5)*x(2) + 51;
            else
                % 2ая функция
                if x(1) <= 0 || x(2) <= 0
                    y = 500;
                else
                    y = x(2)*x(2)*x(2) + 2 * x(1) * x(2) + 1 / (x(1)^2 * x(2)^2);
                end
            end
        end
        
        %% градиент функции
        function y = grad(this, x)
            if this.function_numb == 1
                y = [4 * x(2) + 14 * x(1) + 6 * sqrt(5), 4 * x(1) + 8 * x(2) - 12 * sqrt(5)];
            else
                y = [2 * x(2) - 1 / (2 * sqrt(x(1)^3) * sqrt(x(2))), 3 * x(2)^2 + 2 * x(1) - 1 / (2 * sqrt(x(2)^3) * sqrt(x(1))) ];
            end           
        end
        
        %%
        function y = psi_k(this, hi)
            y = this.f(this.phi_xk + hi * this.phi_pk);
        end
        
        %% сброс счетчика вызова функции
        function reset(this)
            this.count = 0;
        end
                
        %% построение графика функции в виде семейства линий уровня минимизируемых функций
        function draw(this, n)

            if this.function_numb == 1
                a1 = -6;
                b1 = 2;
                a2 = -1;
                b2 = 10;
            else
                a1 = -1;
                b1 = 15;
                a2 = -1;
                b2 = 4;
            end
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
            
            %figure(2);
            %mesh(x2,x1,f);
            
            figure(1);
            levels = -30:5:50;
            [C,h] = contour(x2, x1, f, levels);
            
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
        
                     
        %% метод сопряженных градиентов
        function [resx, resf] = conjugate_gradient(this, x0, eps)
            fprintf('==== Метод минимизации сопряженных градиентов ====\nБазовая точка: x0=[%f, %f]\nТочность: eps=%f;\n\n', [x0(1),x0(2),eps]);
            this.reset();
            w = [];
            xk = x0;
            this.draw_point(xk(2), xk(1));
            fprintf('[%f, %f]\n', [xk(1),xk(2)]);
            k = 1;
            w(k, :) = -this.grad(xk);
            
            if abs(norm(w(1, :))) > eps
                pk = w(1, :);
                while true
                    this.phi_xk = xk;
                    this.phi_pk = pk;
                    hi_k = fminsearch(@this.psi_k, 0, optimset('TolX', this.psi_eps));
                    
                    xk_t = xk + hi_k * pk;
                    
                    w(k+1, :) = - this.grad(xk_t);
                    
                    if abs(norm(w(k+1, :))) < eps
                        this.draw_point(xk_t(2), xk_t(1));
                        this.draw_line(xk, xk_t);
                        fprintf('[%f, %f]\n', [xk_t(1),xk_t(2)]);
                        resx = xk_t;
                        resf = this.f(resx);
                        fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n\n', [resx(1),resx(2),resf,this.count]);
                        return;
                    end
                    
                    this.draw_line(xk, xk_t);
                    
                    check_f = abs(this.f(xk) - this.f(xk_t)) < eps;
                    if check_f
                        resx = xk_t;
                        resf = this.f(xk_t);
                        fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n\n', [resx(1),resx(2),resf,this.count]);
                        return;
                    end
                    
                    xk = xk_t;
                    this.draw_point(xk(2), xk(1));
                    
                    fprintf('[%f, %f]\n', [xk(1),xk(2)]);
                    if mod(k, this.conjugate_reset) == 0
                        pk = w(k+1, :);
                    else
                        gamma_k = norm(w(k+1,:))^2 / norm(w(k,:))^2;
                        pk = gamma_k * pk + w(k+1,:);
                    end
                    k = k+1;
                end
            end
        end
        
        %% матрица Гессе
        function H = hesse_matrix(this, x)
           h = 1e-2;
           t = this.f(x);
           tt = (this.f(x + h*[1 1]) - this.f(x + h*[1 -1]) - this.f(x + h*[-1 1]) + this.f(x + h*[-1 -1])) / (4 * h^2);
           H = [(this.f(x - h*[1 0]) - 2*t + this.f(x + h*[1 0]))/ h^2, ...
                tt; ...
                tt, ...
                (this.f(x - h*[0 1]) - 2*t + this.f(x + h*[0 1]))/ h^2];
        end
        
        %% метод Ньютона с конечно-разностной аппроксимацией производных
        function [resx, resf] = newton_with_approximation(this, x0, eps)
            fprintf('==== Метод Ньютона с конечно-разностной аппроксимацией производных ====\nБазовая точка: x0=[%f, %f]\nТочность: eps=%f;\n\n\n', [x0(1),x0(2),eps]);
            this.reset();
            xk = x0;
            while true
                this.draw_point(xk(2), xk(1));
                fprintf('[%f, %f]\n', [xk(1),xk(2)]);
                g = this.grad(xk);
                hesse = this.hesse_matrix(xk);
                check_hesse = hesse(1, 1) > 0 && det(hesse) > 0;
                
                
                if abs(norm(g)) < eps
                    resx = xk;
                    resf = this.f(xk);
                    
                    if ~check_hesse
                        disp('Матрица Гессе не положительно определена! Необходимо провести дополнительный анализ в окрестности x*');
                    end
                    
                    fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n\n', [resx(1),resx(2),resf,this.count]);
                    return;
                else
                    H = hesse;
                    if ~check_hesse
                        disp('Матрица Гессе не положительно определена! Необходимо ее изменить!');
                    end
                    t = H \ (-g');
                    xkt = xk + t';
                    this.alpha = 1;
                    while this.f(xkt) > this.f(xk)
                        this.alpha = this.alpha / 2;
                        xkt = xk + this.alpha * t';                        
                    end
                    this.draw_line(xk, xkt);
                    if abs(mean(xk - xkt)) < eps
                        resx = xkt;
                        resf = this.f(xkt);
                        fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n\n', [resx(1),resx(2),resf,this.count]);
                        return;
                    end
                    
                    dif = abs(this.f(xk) - this.f(xkt));
                    if dif < eps
                        resx = xkt;
                        resf = this.f(xkt);
                        fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n\n', [resx(1),resx(2),resf,this.count]);
                        return;
                    end
                    
                    xk = xkt;
                end
            end
        end
        
        %% ДФП
        function [resx, resf] = DFP(this, x0, eps)
            fprintf('==== Метод Девидсона-Флетчера-Пауэла (ДФП) ====\nБазовая точка: x0=[%f, %f]\nТочность: eps=%f;\n\n\n', [x0(1),x0(2),eps]);
            this.reset();
            
            w = [];
            xk = x0;
            A = eye(this.n);
            k = 1;
            
            while true
                this.draw_point(xk(2), xk(1));
                w(1, :) = -this.grad(xk);
                
                
                if abs(norm(w(k, :), 2)) < eps
                    resx = xk;
                    resf = this.f(xk);
                    fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n\n', [resx(1),resx(2),resf,this.count]);
                    return;
                end
                
                pk = (A*w(k,:)')';
                this.phi_xk = xk;
                this.phi_pk = pk;
                hi_k = fminsearch(@this.psi_k, 0, optimset('TolX', this.psi_eps));
                
                t = xk + hi_k * pk;
                
                w(k+1,:) = - this.grad(t);
                
                if mod(k, this.conjugate_reset) == 0
                    A = eye(2);
                else
                    dxK = (t-xk)';
                    this.draw_line(xk, t);
                    if abs(this.f(xk) - this.f(t)) < eps
                        resx = t;
                        resf = this.f(t);
                        fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n\n', [resx(1),resx(2),resf,this.count]);
                        return;
                    end
                    
                    %dif = abs(this.f(xk) - this.f(t));
                    %if dif < eps
                    %    resx = t;
                    %    resf = this.f(t);
                    %    fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n\n', [resx(1),resx(2),resf,this.count]);
                    %    return;
                    %end
                    xk = t;
                    
                    dwK = (w(k+1, :) - w(k,:))';
                    
                    A = A - dxK*(dxK')/sum(dwK.*dxK) - A * dwK * (dwK') * (A') / sum((A*dwK).*dwK);
                    
                end
                k = k + 1;
                
            end
        end
        
        %% использование возможностей Optimization Toolbox Matlab
        function [resx, resf] = fminsearch(this, x0, eps)
            this.reset();
            fprintf('==== Optimization Toolbox Matlab ====\nБазовая точка: x0=[%f, %f]\nТочность: eps=%f;\n\n', [x0(1),x0(2),eps]);
            [resx, resf] = fminsearch(@this.f, x0, optimset('TolX', eps));
            fprintf('>>>> Достигнута требуемая точность.\nТочка минимума: xmin=[%f, %f], fmin=%f;\nКол-во вычислений функции: count=%d.\n\n', [resx(1),resx(2),resf, this.count]);
        end
    end
end