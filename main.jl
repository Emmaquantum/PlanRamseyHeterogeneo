# ----------------------
# Proyecto Final
# Análisis Númerico I
# Universidad Nacional de Colombia
# 
# Manuel Armando Alvarado Pinzón
# Julian Sierra Salamanca
# ----------------------

# ----------------------
# Librerías
# ----------------------
using LinearAlgebra
using Interpolations
using Plots

# ----------------------
# Parámetros del modelo
# ----------------------
const α = 0.36         # Participación del capital
const β = 0.995        # Descuento
const δ = 0.005        # Depreciación
const σ = 2.0          # CRRA
const τ = 0.02         # Impuesto inicial
const ρ = 0.25         # Tasa de reemplazo desempleo

# ----------------------
# Tolerancias
# ----------------------
const tol₁   = 1e-7    # Tolerancia Golden-Section
const tolᵥ   = 1e-3    # Tolerancia Bellman
const negval = -1e10   # Penalización consumo negativo
const ε      = 0.05    # Step para bracket boundary
const nᵢₜ    = 100     # Iteraciones máximas Bellman

# ----------------------
# Grilla
# ----------------------
const aₘᵢₙ = -2.0
const aₘₐₓ = 3000.0
const nₐ   = 201
const a    = range(aₘᵢₙ, stop = aₘₐₓ, length = nₐ)

# ----------------------
# Parámetros del Golden-Section
# ----------------------
const r₁ = 0.61803399
const r₂ = 1 - r₁

# ----------------------
# Matriz de transición de empleo
# ----------------------
const prob = [0.9565 0.0435; 0.5 0.5]

# ----------------------
# Matriz de transición de empleo
# ----------------------
include("utility.jl")
include("bellman.jl")
include("golden_section.jl")

# ------------------------------------------------------------------
# 1) Resolver eigen(prob') para autovectores con valor propio de 1
# ------------------------------------------------------------------
eig = eigen(prob')
idx = findall(x -> abs(x - 1) < 1e-8, eig.values)[1]
p₁  = eig.vectors[:, idx]
p₁ .= p₁ ./ sum(p₁)

# Masa de empleados n₀ = p₁[1]
n₀ = p₁[1]

# ------------------------------------------------------------------
# 2) Calcular k₀, w₀, r₀ (q=0)
# ------------------------------------------------------------------
k₀ = (α / (1/β - 1 + δ))^(1/(1 - α)) * n₀
w₀ = (1 - α) * (k₀ / n₀)^α
r₀ = α * k₀^(α - 1) * n₀^(1 - α) - δ
b  = ρ * w₀

# ------------------------------------------------------------------
# 3) Grilla de activos
# ------------------------------------------------------------------
r, w = r₀, w₀

# Funciones de valor y políticas
v    = zeros(nₐ, 2)
aopt = zeros(nₐ, 2)
copt = zeros(nₐ, 2)

# ------------------------------------------------------------------
# 4) Valores iniciales
# ------------------------------------------------------------------
for i in 1:nₐ
    yₑ = (1 - τ) * r * a[i] + (1 - τ) * w
    yᵤ = (1 - τ) * r * a[i] + b
    v[i, 1] = (σ == 1 ? log(yₑ) : yₑ^(1 - σ)/(1 - σ)) / (1 - β)
    v[i, 2] = (σ == 1 ? log(yᵤ) : yᵤ^(1 - σ)/(1 - σ)) / (1 - β)
end

# --------------------------------
# 5) Iteración de Bellman
# --------------------------------
crit = Inf
iter = 0
while iter < nᵢₜ && crit > tolᵥ
    v_old = copy(v)
    v_interp = [LinearInterpolation(a, v_old[:, j], extrapolation_bc = Line()) for j in 1:2]
    for e in 1:2, i in 1:nₐ
        a₀ = a[i]
        vbest = negval
        ax, bx, cx = a[1], a[1], a[end]
        for j in 1:nₐ
            c0 = (1 + (1 - τ) * r) * a₀ + (e == 1 ? (1 - τ) * w : b) - a[j]
            if c0 > 0
                val = utility(c0) + β * (prob[e, 1] * v_interp[1](a[j]) + prob[e, 2] * v_interp[2](a[j]))
                if val > vbest
                    vbest = val
                    ax, bx, cx = a[max(j - 1, 1)], a[j], a[min(j + 1, nₐ)]
                end
            end
        end
        f(x) = bellman(a₀, x, e, v_interp)
        if ax == bx
            bxp = ax + ε * (a[2] - a[1])
            aopt[i, e] = bellman(a₀, bxp, e, v_interp) < bellman(a₀, ax, e, v_interp) ? a[1] : GoldenSectionMax(f, ax, bxp, cx, tol₁)
        elseif bx == cx
            bxp = cx - ε * (a[end] - a[end - 1])
            aopt[i, e] = bellman(a₀, bxp, e, v_interp) < bellman(a₀, cx, e, v_interp) ? a[end] : GoldenSectionMax(f, ax, bxp, cx, tol₁)
        else
            aopt[i, e] = GoldenSectionMax(f, ax, bx, cx, tol₁)
        end
        v[i, e] = bellman(a₀, aopt[i, e], e, v_interp)
    end
    crit = maximum(abs.(v .- v_old))
    iter += 1
end

# --------------------------------
# 6) Política de consumo
# --------------------------------
for i in 1:nₐ
    copt[i, 1] = (1 + (1 - τ) * r) * a[i] + (1 - τ) * w - aopt[i, 1]
    copt[i, 2] = (1 + (1 - τ) * r) * a[i] + b - aopt[i, 2]
end

# --------------------------------
# 7) Grafica función de ahorro neto
# --------------------------------
savingsₑ = aopt[:, 1] .- a
savingsᵤ = aopt[:, 2] .- a

plot(a, savingsₑ, label = "Empleado", linewidth = 2)
plot!(a, savingsᵤ, label = "Desempleado", linewidth = 2)
xlims!(aₘᵢₙ, 1100)
xlabel!("a"); ylabel!("a' − a")
title!("Función estacionaria de ahorro")

# --------------------------------
# 8) Grafica distribución estacionaria de empleo
# --------------------------------
states = ["Empleado", "Desempleado"]
bar(states, p₁, legend = false)
ylabel!("Probabilidad")
title!("Distribución estacionaria de estados de empleo")