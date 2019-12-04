clear all;

######################### Метод для определения матрицы коэффициентов сплайнов ##############################
function matrix_coef = build_matrix_of_coeff(x, ksi)
    ## x - сетка задания функции
    ## ksi - сетка задания сплайна
    iter = 0;
        
    for n = 1:length(x)
        for k = 2:length(ksi)
            if (n == 1)
                S(n,k-1) = 0;
                if (k == length(ksi))
                    S(n,k-1) = ((x(n) - ksi(end))^2)/2;
                end
            else
                if (k == length(ksi))
                    S(n,k-1) = ((x(n) - ksi(end))^2)/2;
                else
                    if (x(n) >= ksi(k) && x(n) <= ksi(k+1))
                        iter = k;                        
                    end
                end
            end
        end
    
        
        if (iter ~= 1)
            for k = 2:iter
                S(n,k-1) = ((x(n) - ksi(k-1))^2 - (x(n) - ksi(k))^2)/2;
            end
        end
    end
    vec_ones = ones(1, length(x));
    matrix_coef = [vec_ones' x' S];
end

########################### Метод для формирования матрицы регуляризации #################################
function omega = get_regularization_matrix(x, ksi)
    ## x - сетка задания функции
    ## ksi - сетка задания сплайна
    
    A(1:length(ksi)-1) = 2;
    A(1) = A(end) = 1;
    B(1:length(ksi)-2) = -1;
    diag_matrix = [zeros(length(ksi)-1,1) zeros(length(ksi)-1,1) (diag(A) + diag(B, -1) + diag(B, 1))];
    
    omega = [zeros(1,length(ksi)+1); zeros(1,length(ksi)+1); diag_matrix];
end

################################ Метод для постороения сплайна ###########################################
function Spline = build_spline(function_vec, number_of_nodes, l, alpha)
    ## function_vec    - вектор значений функции
    ## l               - длина интервала, на котором задана функция
    ## number_of_nodes - количество узлов сплайна
    ## alpha           - параметр регуляризации, который зависит от погрешности, с которой измеряется
    ##                 - заданная функция
    
    N = length(function_vec)-1; %количество участков разбиения
    K = number_of_nodes - 1; %количество участков разбиения сплайна
    h = l/N; %расстояние между соседними узлами сетки функции
    k = l/K; %расстояние между соседними узлами сплайна

    x = 0:h:l; %сетка задания функции
    ksi = 0:k:l; %сетка задания сплайна

    S = build_matrix_of_coeff(x, ksi); %матрица коэффициентов для сплайнов 
    
    omega = get_regularization_matrix(x, ksi); %матрица регуляризации

    A = alpha*omega + S'*S;
    F = S'*function_vec';
    param_of_spline = [A\F];
        
    Spline = S*param_of_spline;
end




######################################## задание тестовых данных #########################################
#---------------------------------------------------------------------------------------------------------
x = 0.01:5*pi/100:5*pi;
f = sin(x)./x; # функция, которую необходимо сгладить
K = length(f); # количество узлов сплайна
l = 10; # длина отрезка, на котором определена функция
alpha = 0.5; # параметр регуляризации

for i = 1:length(f)
    f(i) = f(i) + normrnd(0, 1/10); # добавляем к значениям функции случайную величину
end
#---------------------------------------------------------------------------------------------------------

Spline = build_spline(f, K, l, alpha);

plot(x,Spline, '*b', 'LineWidth', 2, x,f, 'or','LineWidth', 2)
grid on;
    
#print -djpg figure2.jpg
