using StaticCompiler
using LinearAlgebra
using Pkg
Pkg.add("SymmetricMatrices")

function find_max_nondiagonal(A::Symmetric{Float64, Matrix{Float64}})
    n::Int = size(A, 1)
    max_val::Float64 = abs(A[1, 2])
    max_i:: Int = 1
    max_j:: Int = 2
    for i in 1:n
        for j in i+1:n
            current_abs = abs(A[i, j])
            if current_abs > max_val
                max_val = current_abs
                max_i, max_j = i, j
                # alpha < beta
            end
        end
    end
    return abs(A[max_i, max_j]), max_i, max_j
end

function Jacobi_method(A::Symmetric{Float64, Matrix{Float64}}, eps::Float64, maxcount::Int)
    @assert size(A, 1) >= 3
    n:: Int = size(A, 1)
    A_next::Matrix{Float64} = zeros(n, n) 
    # Q_transposed вместо Q для кеш-локальности
    Q_transposed::Matrix{Float64} = Matrix{Float64}(I, n, n)
    Column_alpha_elem::Float64 = 0
    Column_beta_elem::Float64 = 0
    alpha:: Int = 0
    beta:: Int = 0
    max_elem:: Float64 = 0
    varphi:: Float64 = 0
    max_elem, alpha, beta =  find_max_nondiagonal(A)
    count:: Int = 0

    println("Итерация $(count), МаксНедиагЭлем: $(max_elem) из A[$(alpha), $(beta)]")
    display(A)
    #display(Q_transposed)
    while max_elem > eps && count < maxcount
        # arcctg(alpha) = arctan(1/alpha)
        varphi = 1/2 * atan((2*A[alpha, beta])/(A[alpha, alpha]-A[beta, beta]))
        #for j in 1:n
        #    Column_alpha_elem = Q_transposed[alpha, j]*cos(varphi)
        #                      + Q_transposed[beta, j] *sin(varphi)
        #    Column_beta_elem  = Q_transposed[alpha, j]*-sin(varphi)
        #                      + Q_transposed[beta, j] * cos(varphi)
        #    Q_transposed[alpha, j] = Column_alpha_elem
        #    Q_transposed[beta, j]  = Column_beta_elem
        #end
        for i in 1:n
            for j in i:n
                if i == alpha
                    if j == alpha
                        A_next[alpha, alpha] = A[alpha, alpha]*cos(varphi)^2
                             + A[beta, beta] * sin(varphi)^2
                             +2*A[alpha, beta]*sin(varphi)*cos(varphi)
                    elseif j == beta
                        A_next[alpha, beta] = 0
                    else
                        A_next[alpha, j] = A[alpha, j]*cos(varphi)+A[beta, j]*sin(varphi)
                    end
                elseif i == beta
                    if j == alpha
                        @assert(false)
                    end
                    # Это невозможный случай от того, что мы работаем с половиной
                    # матрицы, где i < j; alpha < beta
                    #elseif j == beta
                    if j == beta
                        A_next[beta, beta] = A[alpha, alpha]*sin(varphi)^2
                           + A[beta, beta] * cos(varphi)^2
                           -2*A[alpha, beta]*sin(varphi)*cos(varphi)
                    else
                        A_next[beta, j] = -A[alpha, j]*sin(varphi)+A[beta, j]*cos(varphi)
                    end
                else
                    if j == alpha
                        A_next[i, alpha] = A[alpha, j]*cos(varphi)+A[beta, j]*sin(varphi)
                    elseif j == beta
                        A_next[i, beta] = A[alpha, j]*-sin(varphi)+A[beta, j]*cos(varphi)
                    else
                        A_next[i, j] = A[i, j]
                    end
                end
            end
        end

        @assert(abs(A_next[alpha, alpha]*A_next[alpha, alpha]+
                A_next[beta, beta]*A_next[beta, beta] -
                A[alpha, alpha]*A[alpha, alpha] -
                A[beta, beta]*A[beta, beta] -
                2*A[alpha, beta]*A[alpha, beta])<0.00001)
        A = Symmetric(A_next, :U)
        max_elem, alpha, beta =  find_max_nondiagonal(A)
        count += 1
        if count < 3#% 10 == 0
            println("Итерация $(count), МаксНедиагЭлем: $(max_elem), varphi=$(varphi)")
            display(A)
          #  display(Q_transposed)
        end
    end
    println("Число итераций: $(count)")
    return A, Q_transposed, count
end

N::Int = 13
beta::Float64 = 1+0.1*(52-53)
varepsilon::Float64 = 0.01
println("Задание: найти приближённое полное решение спектральной задачи для матрицы A")
println("Останов на том шаге итерации, когда максимальная по модулю внедиагональная компонента (aka элемент не на главной диагонали) преобразованной матрицы станет меньше \u03b5 = $(varepsilon).")
println("Проверить найденные приближённые собственные векторы и отвечающие им собственные значения матрицы A, проверив соответствующие приближённые равенства A*q+i==lambda*q_i для любого i in 1,4 с указанием погрешности (отдельно показать вектор A*q_i, отдельно вектор lambda_i*q_i и, затем, вектор A*q_i - lambda_i*q_i")
A::Matrix{Float64} = Matrix{Float64}([
                                    10*beta       1       2       3;
                                          1 10*beta      -3      -2;
                                          2      -3 10*beta      -1;
                                          3      -2      -1 10*beta;
])
println("A:")
display(A)
S::Symmetric{Float64, Matrix{Float64}} = Symmetric(A, :U)

count::Int = 0
new_A::Symmetric{Float64, Matrix{Float64}} = Symmetric(zeros(size(A,1), size(A,1)), :U)
Q_transposed::Matrix{Float64} = zeros(size(A,1), size(A,1))
new_A, Q_transposed, count  = Jacobi_method(S, varepsilon, 10000)
println("Полученные A и Q_transposed:")
display(new_A)
display(Q_transposed)
println("Проверка равенств A*q = lambda*q")
n::Int = size(A, 1)
#A_q::Vector{Float64} = zeros(n)
#l_q:: Vector{Float64} = zeros(n)
for i in 1:n
    println("$i: Собственное значение: $(new_A[i, i]), Собственный вектор: $(Q_transposed[i, :])")
    A_q = A * Q_transposed[i, :]
    l_q = new_A[i, i] * Q_transposed[i, :]
    println("A*q = $(A_q), lambda*q = $(l_q)")
    println("A*q-lambda*q=$(A_q-l_q)")
end


