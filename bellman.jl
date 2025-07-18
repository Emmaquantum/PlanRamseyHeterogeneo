# ----------------------
# Función de Bellman dado a₀, a₁, estado e y un interpolador de v
# ----------------------
function bellman(a₀, a₁, e, v_interp)
    c = (1 + (1 - τ) * r) * a₀ + (e == 1 ? (1 - τ) * w : b) - a₁
    if c < 0
        return negval
    else
        utility(c) + β * (prob[e, 1] * v_interp[1](a₁) + prob[e, 2] * v_interp[2](a₁))
    end
end