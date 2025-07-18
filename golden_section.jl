# ----------------------
# Búsqueda de máximo con Golden Section
# ----------------------
function GoldenSectionMax(f, ax, bx, cx, tol)
    x0, x3 = ax, cx
    if abs(cx - bx) <= abs(bx - ax)
        x1, x2 = bx, bx + r₂ * (cx - bx)
    else
        x2, x1 = bx, bx - r₂ * (bx - ax)
    end
    f1, f2 = f(x1), f(x2)
    while abs(x3 - x0) > tol * (abs(x1) + abs(x2))
        if f2 > f1
            x0, x1, f1 = x1, x2, f2
            x2 = r₁ * x1 + r₂ * x3
            f2 = f(x2)
        else
            x3, x2, f2 = x2, x1, f1
            x1 = r₁ * x2 + r₂ * x0
            f1 = f(x1)
        end
    end
    f1 > f2 ? x1 : x2
end