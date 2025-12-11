#import Pkg; Pkg.add("ForwardDiff")
using LinearAlgebra
using Plots
using ForwardDiff

function Jacobi_method(A, alpha, beta, n, maxerror = 0.01, maxcount = 10000)
    println("A: ")
    display(A)
    println("alpha: $alpha, beta: $beta, n: $n")
    Matrix_size = size(F)
    @assert beta < n
    @assert alpha < beta
    @assert n == Matrix_size[1]
    @assert n == Matrix_size[2]
    @assert n >= 3
    
    # A[beta][alpha] = max
    # Q = (cos 0 -sin 0)
    #     (0   1    0 0)
    #     (sin 0  cos 0)
    #     (0   0    0 1)
    # varphi: компонента b^alpha_beta матрицы B = (b^i_j) = ^T Q(alpha, beta, varphi)*A[0]*Q(alpha, beta, varphi)=A[1] нулевая
    # b_alpha^alpha = A_a^a cos^2 varphi + A_b^b sin^2 varphi + 2 A_beta^alpha sin varphi cos varphi
    # b_b^b = a_a^a sin^2 f + a^b_b cos^2 f - 2 a^a_b sin f cos f
    # 0 = b^a_b = b^b_a = -(a^a_a - a^b_b) cin f cos f + a^a_b (cos^2 f - sin^2 f)
    # b^a_k = a^a_k cos f + a^b_k sin f
    # b^k_b = -a^k_a sin f + a^k_b cos f
    # b^k_m - a^k_m
    if Matrix_size[1] != Vector_size || Matrix_size[1] != Matrix_size[2] 
        error("Размерности не совпадают")
    end
    if x == 0
        x = ones(Vector_size)
    else
        if Vector_size != size(x, 1) 
            error("Размерности не совпадают")
        end
    end
    println("x0: $x")
    norm_F = norm(F, Inf)
    println("Чебышевская норма матрицы F: $(norm_F)")
    if norm_F >= 1
        error("Чебышевская норма >= 1")
    end

    r0 = F * x + g - x
    norm_r0 = norm(r0, Inf)
    norm_x0 = norm(x, Inf)
    norm_g = norm(g, Inf)

    t_counter = 0
    t_error = (norm_F)/(1-norm_F) * norm_g + norm_F*norm_x0
    println("Верхняя граница погрешности для k=$(t_counter): $(t_error)")
    while true
        t_counter += 1 
        t_error = norm_F * t_error
        println("Верхняя граница погрешности для k=$(t_counter): $(t_error)")
        if t_error <= maxerror || t_counter >= maxcount
            break
        end
    end

    counter = 0
    error_k = norm(F*x + g - x, Inf)
    
    println("Итер: $(counter), погрешность: $(error_k)")
    while error_k > maxerror && counter < maxcount
        x = F * x + g
        counter += 1
        error_k = norm(F*x + g - x, Inf)
        println("Итер: $(counter), погрешность: $(error_k)")
    end
    return x, error_k, counter, t_counter
end

function solve_OLE_via_SI(A, x, b, maxerror)
    println("A:")
    display(A)
    println("b:")
    display(b)
    Matrix_size = size(A)
    Vector_size = size(b, 1) 
    if Matrix_size[1] != Vector_size || Matrix_size[1] != Matrix_size[2] 
        error("Размерности не совпадают")
    end
    
    # F = E - D A
    E = I(Vector_size)
    D = Diagonal(1 ./ diag(A))

    F = E - D * A
    g = D * b
    return simple_iteration_method(F, g, x, maxerror)
end

function Seidel_method(F, g, y_k=0, maxerror = 0.01, maxcount = 10000)
    # СЛАУ: y = F * y + g
    println("F:")
    display(F)
    println("g:")
    display(g)
    Matrix_size = size(F)
    Vector_size = size(g, 1) 
    if Matrix_size[1] != Vector_size || Matrix_size[1] != Matrix_size[2] 
        error("Размерности не совпадают")
    end
    if y_k == 0
        y_k = ones(Vector_size)
    else
        if Vector_size != size(y_k, 1) 
            error("Размерности не совпадают")
        end
    end

    B = LowerTriangular(F) 
    D = Diagonal(diag(F))
    Q = B - D
    P = F - Q
    println("Q:")
    display(Q)
    println("P:")
    display(P)
    norm_F = norm(F, Inf)
    error_k = norm(F*y_k + g - y_k, Inf) 

    counter = 0
    for i in 1:Vector_size  
        print("y_k^$(i) = g^$(i)")
        for j in 1:Vector_size  
            print(" + f^$(i)_$(j) *") 
            (j < i ? print("y_k^$j") : print("y_{k-1}^$j"))
        end
        println("")
    end

    y_k1 = zeros(Vector_size)
    while error_k > maxerror && counter < maxcount
        for i in 1:Vector_size  
            y_k1[i] = g[i]
            for j in 1:Vector_size  
               y_k1[i] += F[i,j] * (j < i ? y_k1[j] : y_k[j])
            end
        end
        y_k = copy(y_k1)  
        counter += 1
    # СЛАУ: y = F * y + g
        error_k = norm(F*y_k + g - y_k, Inf)  
        println("Итер: $(counter), оценка погрешности: $(error_k)")
    end
    return y_k, error_k, counter
end

function solve_OLE_via_Seidel(A, x, b, maxerror)
    Matrix_size = size(A)
    Vector_size = size(b, 1) 
    if Matrix_size[1] != Vector_size || Matrix_size[1] != Matrix_size[2] 
        error("Размерности не совпадают")
    end
    
    # F = E - D A
    E = I(Vector_size)
    D = Diagonal(1 ./ diag(A))

    F = E - D * A
    g = D * b
    return Seidel_method(F, g, x, maxerror)

end

N = 13
n = 53
beta = 1-0.03*(50-n)
println("beta: $(beta)")

A = [10*beta       1       2       3;
           1 10*beta       3      -2;
           2      -3 10*beta       1;
           3       2      -1  10*beta]

x = ones(4)
b = A*x
x = zeros(4)  
x, error_k, counter, t_counter = solve_OLE_via_SI(A, x, b, 0.01)
println("x: $(x)")
println("Погрешность: $(error_k), Число итераций: $(counter), Теоретическая оценка максимального числа итераций: $(t_counter)")
x = zeros(4)
x, error_k, counter = solve_OLE_via_Seidel(A, b, x, 0.01)
println("Погрешность: $(error_k), Число итераций: $(counter)")

function Newton_method(f, a, b, maxerror, maxcount = 1000)
    # x = x - f(x)/f'(x)      
    println("f: $f, a: $a, b: $b, maxerror: $maxerror") 
    x = (a+b)/2

    counter = 0
    error_k = abs(f(x))
    println("Итер: $(counter), погрешность: $(error_k)")
    while error_k > maxerror && counter < maxcount
        x = x - f(x)/(ForwardDiff.derivative(f, x)) 
        counter += 1
        error_k = abs(f(x))
        println("Итер: $(counter), x: $x, погрешность: $(error_k)")
    end
    return x, error_k, counter
end

function sekuschaja(f, a, b)
    return (a-b)*f(a)/(f(a)-f(b))
end


# x0 = a
function working_formula1(f, x, b)
    return (b-x)*f(x)/(f(b)-f(x))
end

# x0 = b
function working_formula2(f, x, a)
    return (x-a)*f(x)/(f(x)-f(a))
end

d2f = t -> ForwardDiff.derivative(y -> ForwardDiff.derivative(f, y), t)
d3f = t -> ForwardDiff.derivative(d2f, t)
function method_sek(f, a, b, maxerror = 0.01, maxcount = 1000)
    # x = x - f(x0/varphi(x)
    # vvarphi(x) = (x - x_k-1)f(x)/(f(x)-f(x_k-1)
    x = (a+b)/2
    counter = 0
    error_k = abs(f(x))
    println("f(a)*f(b): $(f(a)*f(b))")
    can_work1 = false
    can_work2 = false
    if (f(a)*f(b) < 0)
        s_der = d2f((a+b)/2)
        t_der = d3f((a+b)/2)
        s_der0 = s_der + t_der* -((a+b)/2)
        s_der1 = s_der + t_der* ((a+b)/2)
        the_sign = sign(s_der)
        can_work = the_sign == sign(s_der0) && the_sign == sign(s_der1)
        if can_work
            println("f'' считается одного знака, знак: $the_sign")
            println("f(b)*f'': $(f(b)*s_der)")
            if f(b)*s_der > 0
                println("Целесообразно использовать метод, использующий рабочую формулу 1")
                can_work1 = true
            else
                println("Использовать метод, использующий рабочую формулу 1, нецелесообразно")
                can_work1 = false
            end

            if f(a)*s_der > 0
                println("Целесообразно использовать метод, использующий рабочую формулу 2")
                can_work2 = true
            else
                println("Использовать метод, использующий рабочую формулу 2, нецелесообразно")
                can_work2 = false
            end
        else
            println("f'' принимет разные знаки")
        end
    end

    if can_work1
        x = a
        while error_k > maxerror && counter < maxcount
            x = x - (working_formula1(f, x, b)) 
            counter += 1
            error_k = abs(f(x))
            println("Итер: $(counter), x: $x, погрешность: $(error_k)")
        end
        return x, error_k, counter
    end
    
    if can_work2
        x = b
        while error_k > maxerror && counter < maxcount
            x = x - (working_formula2(f, x, a)) 
            counter += 1
            error_k = abs(f(x))
            println("Итер: $(counter), x: $x, погрешность: $(error_k)")
        end
        return x, error_k, counter
    end
    x_prev = (a+b)/2
    println("Двушаговый метод секущих")
    x = x_prev
    x = x - f(x)/(ForwardDiff.derivative(f, x)) 
    counter += 1
    error_k = abs(f(x))
    while error_k > maxerror && counter < maxcount
        x_new = x - f(x)* (x-x_prev)/(f(x) - f(x_prev)) 
        println("x_new: $x_new, x: $x, x_prev: $x_prev")
        x_prev = x
        x = x_new
        counter += 1
        error_k = abs(f(x))
        println("Итер: $(counter), погрешность: $(error_k)")
    end
    return x, error_k, counter
end

N = 13
n = 53
alpha = 0.003*(n-50)
f(t)  =  (N + 5.2 + (-1)^N * alpha) * t^3 -
           (2*N^2 + 10.4*N + (-1)^(N+1)*alpha) * t^2 -
           N^2*(N + 5.2)*(t-2*N) + (-1)^N*alpha
#f(x) =  let N=N, alpha=alpha
#    x -> (N + 5.2 + (-1)^N * alpha) * x^3 -
#         (2*N^2 + 10.4*N + (-1)^(N+1)*alpha) * x^2 -
#         N^2*(N + 5.2)*(x-2*N) + (-1)^N*alpha
#end
#f(x) = (N + 5.2 + (-1)^N * alpha) * x^3
#            - (2*N^2 + 10.4*N + (-1)^(N+1)*alpha) * x^2 
#            - N^2*(N + 5.2)*(x-2*N)
#            + (-1)^N*alpha

#plot(t -> f(t), -50, 50)
#display(plot)
#savefig("График.png")

#plot(t -> f(t), 13, 28)
#display(plot)
#savefig("График2.png")

#plot(t -> f(t), -20, -10)
#display(plot)
#savefig("График3.png")



x, error_k, counter = Newton_method(f, -20, -10, 0.0001)
println("Метод Ньютона (касательных):\nx: $x, погрешность: $error_k, Число итераций: $counter")

x, error_k, counter = method_sek(f, 24, 30, 0.0001)
println("Метод секущих:\nx: $x, погрешность: $error_k, Число итераций: $counter")

function divide_et_vince(f, a, b, maxerror, maxcount = 1000)
    # x = x - f(x)/f'(x)      
    println("f: $f, a: $a, b: $b, maxerror: $maxerror") 
    x = (a+b)/2
    @assert f(a)*f(b) < 0 "Если f(a)*f(b)>=0, то метод не применим"
    counter = 0
    error_k = abs(f(x))
    println("Итер: $(counter), погрешность: $(error_k)")
    upper = b
    bottom = a
    the_sign = sign(f(bottom))
    while error_k > maxerror && counter < maxcount
        if sign(f(x)) != the_sign
            upper = x
        else
            bottom = x
        end
        x = (upper + bottom)/2
        counter += 1
        error_k = abs(f(x))
        println("Итер: $(counter), погрешность: $(error_k)")
        if counter % 10 == 0
            @assert sign(f(upper)) != sign(f(bottom)) "Если f(a)*f(b)>=0, то метод не применим"
        end
    end
    return x, error_k, counter
end


x, error_k, counter = divide_et_vince(f, 0, 21, 0.0001)
println("Метод средних:\nx: $x, погрешность: $error_k, Число итераций: $counter")
#function method_absciss(f, a, b, x=0, maxerror = 0.01, maxcount = 10000)
#function Newton_method(f, a, b, x=0, maxerror = 0.01, maxcount = 10000)
