using LinearAlgebra

function simple_iteration_method(F, g, x=0, maxerror = 0.01, maxcount = 10000)
    Matrix_size = size(F)
    Vector_size = size(g, 1)
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
    norm_g = norm(g, Inf)

    if norm_F >= 1
        error("Чебышевская норма >= 1")
    end
    
    counter = 0
    error_k = Inf
    t_counter = 1
    norm_Fk = norm_F
    norm_x = norm(x, Inf)
    t_error = error_k
    while t_error > maxerror && t_counter < maxcount
        t_error = norm_Fk/(1-norm_F) * norm_g + norm_F * norm_x
        norm_Fk *= norm_F
        t_counter += 1
        println("Iter: $t_counter, theoretical error: $t_error")
    end

    while error_k > maxerror && counter < maxcount
        error_k = norm(F*x + g - x, Inf) 
        x = F * x + g
        counter += 1
        println("Iter: $(counter), error: $(error_k)")
    end
   # println("error_k: $(error_k), maxerror: $maxerror")
    #println("counter: $counter, maxcount: $maxcount")
    error_k = norm(F * x + g - x, Inf)  
    return x, error_k, counter, t_counter
end

function solve_OLE_via_SI(A, x, b, maxerror)
    Matrix_size = size(A)
    Vector_size = size(b, 1) 
    if Matrix_size[1] != Vector_size || Matrix_size[1] != Matrix_size[2] 
        error("Размерности не совпадают")
    end
    
    # F = E - D A
    E = Matrix(I, Vector_size, Vector_size) 
    D = Diagonal(1 ./ diag(A))

    F = E - D * A
    g = D * b
    return simple_iteration_method(F, g, x, maxerror)
end

function Seidel_method(F, g, y_k=0, maxerror = 0.01, maxcount = 10000)
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
    
    error_k = norm(F*y_k + g - y_k, Inf) 

    y_k1 = zeros(Vector_size)
    counter = 0
    while error_k > maxerror && counter < maxcount
        for i in 1:Vector_size  
            y_k1[i] = g[i]
            for j in 1:Vector_size  
                if j != i
                    y_k1[i] += F[i,j] * (j < i ? y_k1[j] : y_k[j])
                end
            end
        end
        y_k = copy(y_k1)  
        counter += 1
        error_k = norm(F*y_k + g - y_k, Inf)  
        if counter % 10 == 0
            println("Итерация: $(counter), оценка погрешности: $(error_k)")
        end
    end
    return y_k, error_k, counter
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
println("Error: $(error_k), Counter: $(counter), Theoretical counter: $(t_counter)")
