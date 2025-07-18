# ----------------------
# Función de utilidad CRRA/Log
# ----------------------
function utility(c)
    σ == 1 ? log(c) : c^(1 - σ) / (1 - σ)
end