using LinearAlgebra

function cond_manual(A)
    invA = inv(A)
    return opnorm(A, Inf) * opnorm(invA, Inf)
end

function relative_solution_error(A, b, Db)
    c = cond_manual(A)
    return c * (norm(Db) / norm(b)), c
end

function Task1()
    println("Task 1")
    N = 13
    α = 0.03

    A = [
        200*(1+0.5N+α)   200*(1+0.5N)         200*(1+0.5N)
        200.1*(1+0.5N)  199.9*(1+0.5N+α)    200*(1+0.5N)
        199.9*(1+0.5N)  200*(1+0.5N)        200.1*(1+0.5N+α)
    ]

    bval = 200*(3 + 1.5N + α)
    b  = fill(bval, 3)
    Db = [0.01bval, -0.01bval, 0.01bval]

    rel_err, c = relative_solution_error(A, b, Db)

    println("cond(A) = $c")
    println(c > 100 ? "SLE is ill-conditioned" : "SLE is well-conditioned")
    println("Relative error ≤ $rel_err")
    println("───────────────")
end

function Task2()
    println("Task 2")

    N = 13
    α = 0.02
    λ = 0.4 - α
    F(x) = atan(x)

    a = 0
    b = 1
    GRID = 10
    h = (b - a) / GRID

    s = [a + h/2 + i*h for i in 0:GRID-1]

    A = [F(s[i]*s[j]) * (b - a) / GRID for i in 1:GRID, j in 1:GRID]

    Aλ = I + λ * A

    x = ones(GRID)

    bvec = Aλ * x

    imat_cond = cond_manual(Aλ)
    println("cond(E + λA) = $imat_cond")
    println(imat_cond > 100 ? "SLE is ill-conditioned" : "SLE is well-conditioned")

    Db = [(-1)^(i-1) * 0.01 * bvec[i] for i in 1:GRID]

    # оценка относительной погрешности
    rel_err1 = imat_cond * norm(Db) / norm(bvec)
    println("Relative error ≤ $rel_err1")

    scale = bvec .+ Db

    Aλ2 = [Aλ[i,j] / scale[i] for i in 1:GRID, j in 1:GRID]
    b2  = bvec ./ scale
    Db2 = Db    ./ scale

    imat_cond2 = cond_manual(Aλ2)
    println("cond(scaled system) = $imat_cond2")

    rel_err2 = imat_cond2 * norm(Db2) / norm(b2)
    println("Relative error (scaled) ≤ $rel_err2")
    println("───────────────")
end

Task1()
Task2()
