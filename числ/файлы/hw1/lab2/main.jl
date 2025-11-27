import Pkg; Pkg.add("Plots"); Pkg.add("Distributions")

using LinearAlgebra
using Statistics
using Plots
using Distributions

y_lin = [
    17.89, 16.79, 17.60, 20.88, 16.67,
    17.53, 19.27, 16.51, 14.95, 16.47,
    17.66, 16.08, 18.96, 16.49, 15.37,
    15.27, 19.24, 17.43, 17.75, 15.52
]

Z = [
  1.158574 1.194067 1.745872 1.566271 1.825556 1.942503;
  1.238868 1.913419 1.182653 1.044649 1.304209 1.924039;
  1.564043 1.561357 1.070589 1.778954 1.226447 1.824122;
  1.737266 1.798975 1.952239 1.752281 1.247871 1.547960;
  1.364544 1.031220 1.380596 1.688101 1.987396 1.058504;
  1.535295 1.742973 1.580401 1.063356 1.999237 1.425459;
  1.780725 1.306711 1.972594 1.686270 1.582629 1.767235;
  1.135044 1.139164 1.686178 1.220069 1.034577 1.019745;
  1.246498 1.114597 1.079653 1.333415 1.054445 1.156743;
  1.416456 1.349223 1.680380 1.003235 1.471908 1.095523;
  1.611866 1.972991 1.443953 1.014008 1.916990 1.182531;
  1.520585 1.427992 1.464156 1.011505 1.108341 1.981536;
  1.229896 1.304392 1.852107 1.705496 1.725639 1.214820;
  1.726829 1.866756 1.074984 1.098880 1.983154 1.256935;
  1.772790 1.363353 1.227454 1.076754 1.656758 1.675253;
  1.418256 1.072481 1.123447 1.438917 1.059481 1.080325;
  1.119724 1.947356 1.372631 1.635578 1.940580 1.112827;
  1.728446 1.802332 1.365001 1.184759 1.119633 1.880032;
  1.161107 1.359294 1.956206 1.143406 1.491440 1.688437;
  1.963561 1.271859 1.250008 1.193670 1.466262 1.624409
]

m = length(y_lin)
if size(Z,1) != m
    error("Размерности Z и y не совпадают")
end

A_lin = hcat(ones(m), Z) 
p_lin = size(A_lin,2)   

beta_lin = (A_lin' * A_lin) \ (A_lin' * y_lin) 
yhat_lin = A_lin * beta_lin
res_lin = y_lin .- yhat_lin
sigma2_lin = sum(res_lin .^ 2) / (m - p_lin)  
sigma_lin = sqrt(sigma2_lin)

println("=== Задание 1-2.1 (линейная модель) ===")
println("Оценки параметров (x0..x6):")
for (i, b) in enumerate(beta_lin)
    println("  x$(i-1) = $(round(b, digits=6))")
end
println("Оценка дисперсии sigma^2_hat = $(round(sigma2_lin, digits=6)), sigma_hat = $(round(sigma_lin, digits=6))")
println()

plot1 = scatter(1:m, y_lin, label="y (набл.)", legend=:topleft, title="Линейная модель: y и y_hat (N=13)",
    xlabel="i", ylabel="y", size=(1280,720))
plot!(1:m, yhat_lin, label="y_hat (прогноз)", lw=2)
plot2 = scatter(1:m, res_lin, label="остатки", title="Линейная модель: остатки (N=13)", xlabel="i", ylabel="остаток", size=(1280,720))
display(plot1)
display(plot2)


println("=== Редукция линейной модели по значимости (α = 0.05) ===")

function reduce_model(y, Z; alpha = 0.05)
    m, k = size(Z)         
    A = hcat(ones(m), Z)  
    p = size(A,2)        
    
    beta = (A' * A) \ (A' * y)
    yhat = A * beta
    res = y - yhat
    sigma2 = sum(res.^2) / (m - p)
    sigma = sqrt(sigma2)
    
    cov_beta = sigma2 * inv(A' * A)
    
    stderr = sqrt.(diag(cov_beta))
    
    tvals = beta ./ stderr
    
    tcrit = quantile(TDist(m - p), 1 - alpha/2)
    
    significant = abs.(tvals) .> tcrit
    
    println("\nПолная модель: p = $p параметров")
    println("Критическое значение t = $(round(tcrit, digits=4))")
    println("Коэффициенты, t-статистики и значимость:")

    for i in 1:p
        println("x$(i-1): β=$(round(beta[i], digits=5)),  t=$(round(tvals[i], digits=4)),  significant=$(significant[i])")
    end
    
    if all(significant)
        println("\nВсе коэффициенты значимы. Редукция модели невозможна.")
        return beta, sigma, collect(0:k)
    end
    
    idx = findall(significant) 
    A_red = A[:, idx]
    p_red = size(A_red, 2)
    
    println("\nРедуцированная модель: оставлены параметры x^{", join(idx .- 1, ", "), "}")
    
    beta_red = (A_red' * A_red) \ (A_red' * y)
    yhat_red = A_red * beta_red
    res_red = y - yhat_red
    sigma2_red = sum(res_red.^2) / (m - p_red)
    sigma_red = sqrt(sigma2_red)
    
    println("\nИтог после редукции:")
    println("Число параметров: $p_red")
    println("sigma_hat = $(round(sigma_red, digits=5))")
    
    return beta_red, sigma_red, idx .- 1 
end


beta_red, sigma_red, significant_idx = reduce_model(y_lin, Z)
println("\nОставшиеся значимые параметры: x^($(join(significant_idx, ", ")))")


cols = []

for idx in significant_idx
    if idx == 0
        push!(cols, ones(m))
    else
        push!(cols, Z[:, idx])
    end
end

A_red = hcat(cols...)

yhat_red = A_red * beta_red
res_red = y_lin .- yhat_red

plot_red1 = scatter(
    1:m, y_lin,
    label="y (набл.)",
    title="Редуцированная модель: y и y_hat",
    xlabel="i",
    ylabel="y",
    size=(1280, 720)
)

plot!(1:m, yhat_red, label="y_hat (редуц.)", lw=2)
display(plot_red1)
savefig(plot_red1, "plot5.png")

plot_red2 = scatter(
    1:m, res_red,
    label="остатки (редуц.)",
    title="Редуцированная модель: остатки",
    xlabel="i",
    ylabel="остаток",
    size=(1280, 720)
)
display(plot_red2)
savefig(plot_red2, "plot6.png")

y_poly = [
   -4.94, -4.73, -4.53, -4.34, -4.16,
   -4.00, -3.83, -3.66, -3.49, -3.31,
   -3.11, -2.90, -2.67, -2.42, -2.13,
   -1.82, -1.47, -1.09, -0.66, -0.18
]

t = collect(0.0:0.05:0.95) 
if length(t) != length(y_poly)
    error("t и y_poly должны быть одинаковой длины")
end

k_poly = 4
A_poly = hcat([t .^ j for j in 0:k_poly]...)

p_poly = size(A_poly,2) 

beta_poly = (A_poly' * A_poly) \ (A_poly' * y_poly)
yhat_poly = A_poly * beta_poly
res_poly = y_poly .- yhat_poly
sigma2_poly = sum(res_poly .^ 2) / (m - p_poly)
sigma_poly = sqrt(sigma2_poly)

println("=== Задание 1-2.2 (полином 4-й степени) ===")
println("Оценки параметров (x0..x4):")
for (i, b) in enumerate(beta_poly)
    println("  x$(i-1) = $(round(b, digits=6))")
end
println("Оценка дисперсии sigma^2_hat = $(round(sigma2_poly, digits=6)), sigma_hat = $(round(sigma_poly, digits=6))")
println()

tt = range(minimum(t), stop=maximum(t), length=200)
A_tt = hcat([tt .^ j for j in 0:k_poly]...)
yfit_smooth = A_tt * beta_poly

plot3 = scatter(t, y_poly, label="y (набл.)", title="Полином степени 4: y и аппроксимация (N=13)",
    xlabel="t", ylabel="y", legend=:topright, size=(1280,720))
plot!(tt, yfit_smooth, label="аппрокс.", lw=2)
plot4 = scatter(1:m, res_poly, label="остатки", title="Полином степени 4: остатки (N=13)", xlabel="i", ylabel="остаток", size=(1280,720))
display(plot3)
display(plot4)
savefig(plot1, "plot1.png")
savefig(plot2, "plot2.png")
savefig(plot3, "plot3.png")
savefig(plot4, "plot4.png")
