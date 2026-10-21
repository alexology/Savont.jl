# Exact ports of the two statistical routines used by upstream marker inference:
# statrs 0.16.1 (`Binomial::cdf` through the regularized incomplete beta function)
# and fishers_exact 1.0.1 (Langsrud's two-sided Fisher exact test). The operation
# order follows the Rust sources so that p-values agree to the last few ulps.

const _GAMMA_R = 10.900511
const _GAMMA_DK = (
    2.48574089138753565546e-5,
    1.05142378581721974210,
    -3.45687097222016235469,
    4.51227709466894823700,
    -2.98285225323576655721,
    1.05639711577126713077,
    -1.95428773191645869583e-1,
    1.70970543404441224307e-2,
    -5.71926117404305781283e-4,
    4.63399473359905636708e-6,
    -2.71994908488607703910e-9,
)
const _LN_PI = 1.1447298858494001741434273513530587116472948129153
const _LN_2_SQRT_E_OVER_PI = 0.6207822376352452223455184457816472122518527279025978
const _F64_PREC = 0.00000000000000011102230246251565

function _ln_gamma(x::Float64)
    if x < 0.5
        s = _GAMMA_DK[1]
        for t in 1:10
            s += _GAMMA_DK[t + 1] / (t - x)
        end
        return _LN_PI - log(sin(pi * x)) - log(s) - _LN_2_SQRT_E_OVER_PI -
               (0.5 - x) * log((0.5 - x + _GAMMA_R) / ℯ)
    end
    s = _GAMMA_DK[1]
    for t in 1:10
        s += _GAMMA_DK[t + 1] / (x + t - 1.0)
    end
    return log(s) + _LN_2_SQRT_E_OVER_PI + (x - 0.5) * log((x - 0.5 + _GAMMA_R) / ℯ)
end

# approx 0.5 `ulps_eq!` with its defaults (epsilon = eps(Float64), max_ulps = 4).
function _ulps_eq(a::Float64, b::Float64)
    abs(a - b) <= eps(Float64) && return true
    signbit(a) == signbit(b) || return false
    ia = reinterpret(Int64, a); ib = reinterpret(Int64, b)
    return abs(ia - ib) <= 4
end

function _beta_reg(a::Float64, b::Float64, x::Float64)
    a > 0.0 || throw(ArgumentError("a must be positive"))
    b > 0.0 || throw(ArgumentError("b must be positive"))
    0.0 <= x <= 1.0 || throw(ArgumentError("x must be in [0, 1]"))
    bt = x == 0.0 || _ulps_eq(x, 1.0) ? 0.0 :
         exp(_ln_gamma(a + b) - _ln_gamma(a) - _ln_gamma(b) + a * log(x) + b * log(1.0 - x))
    symm_transform = x >= (a + 1.0) / (a + b + 2.0)
    fpmin = floatmin(Float64) / _F64_PREC
    if symm_transform
        a, b, x = b, a, 1.0 - x
    end
    qab = a + b; qap = a + 1.0; qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    abs(d) < fpmin && (d = fpmin)
    d = 1.0 / d
    h = d
    for mi in 1:140
        m = Float64(mi); m2 = m * 2.0
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        abs(d) < fpmin && (d = fpmin)
        c = 1.0 + aa / c
        abs(c) < fpmin && (c = fpmin)
        d = 1.0 / d
        h = h * d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        abs(d) < fpmin && (d = fpmin)
        c = 1.0 + aa / c
        abs(c) < fpmin && (c = fpmin)
        d = 1.0 / d
        del = d * c
        h *= del
        if abs(del - 1.0) <= _F64_PREC
            return symm_transform ? 1.0 - bt * h / a : bt * h / a
        end
    end
    return symm_transform ? 1.0 - bt * h / a : bt * h / a
end

"""
    binomial_test(n, k, p) -> Float64

Upper binomial tail used by upstream `utils::binomial_test`: `1 - P(X <= k)` for
`X ~ Binomial(n, p)`, i.e. `P(X > k)` (the upstream comment says "k or more", but
the computed quantity excludes `k`).
"""
function binomial_test(n::Integer, k::Integer, p::Real)
    n >= 0 && k >= 0 || throw(ArgumentError("n and k must be nonnegative"))
    0 <= p <= 1 || throw(ArgumentError("p must be in [0, 1]"))
    cdf = k >= n ? 1.0 : _beta_reg(Float64(n - k), Float64(k) + 1.0, 1.0 - Float64(p))
    return 1.0 - cdf
end

function _lngamm(zi::Int)
    z = Float64(zi)
    x = 0.0
    x += 0.1659470187408462e-06 / (z + 7.0)
    x += 0.9934937113930748e-05 / (z + 6.0)
    x -= 0.1385710331296526 / (z + 5.0)
    x += 12.50734324009056 / (z + 4.0)
    x -= 176.6150291498386 / (z + 3.0)
    x += 771.3234287757674 / (z + 2.0)
    x -= 1259.139216722289 / (z + 1.0)
    x += 676.5203681218835 / z
    x += 0.9999999999995183
    return log(x) - 5.58106146679532777 - z + (z - 0.5) * log(z + 6.5)
end

_lnfact(n::Int) = n <= 1 ? 0.0 : _lngamm(n + 1)
_lnbico(n::Int, k::Int) = _lnfact(n) - _lnfact(k) - _lnfact(n - k)
_hyper_323(n11, n1_, n_1, n) = exp(_lnbico(n1_, n11) + _lnbico(n - n1_, n_1 - n11) - _lnbico(n, n_1))

mutable struct _HyperState
    n11::Int
    n1_::Int
    n_1::Int
    n::Int
    prob::Float64
    valid::Bool
end

function _hyper0!(s::_HyperState, n11i::Int, n1_i::Int, n_1i::Int, ni::Int)
    if s.valid && (n1_i | n_1i | ni) == 0
        if n11i % 10 != 0
            if n11i == s.n11 + 1
                s.prob *= ((s.n1_ - s.n11) / n11i) * ((s.n_1 - s.n11) / (n11i + s.n - s.n1_ - s.n_1))
                s.n11 = n11i
                return s.prob
            end
            if n11i == s.n11 - 1
                s.prob *= (s.n11 / (s.n1_ - n11i)) * ((s.n11 + s.n - s.n1_ - s.n_1) / (s.n_1 - n11i))
                s.n11 = n11i
                return s.prob
            end
        end
        s.n11 = n11i
    else
        s.n11 = n11i; s.n1_ = n1_i; s.n_1 = n_1i; s.n = ni; s.valid = true
    end
    s.prob = _hyper_323(s.n11, s.n1_, s.n_1, s.n)
    return s.prob
end
_hyper!(s::_HyperState, n11::Int) = _hyper0!(s, n11, 0, 0, 0)

"""
    fishers_exact(table) -> Float64

Two-sided Fisher exact p-value for the 2×2 table `(n11, n12, n21, n22)`, ported
from the `fishers_exact` crate (Langsrud's algorithm, capped at 1).
"""
function fishers_exact(table)
    length(table) == 4 || throw(ArgumentError("table must have four cells"))
    all(x -> 0 <= x <= typemax(Int32), table) ||
        throw(ArgumentError("table cells must be in 0:$(typemax(Int32))"))
    n11, n12, n21, n22 = Int.(Tuple(table))
    n1_ = n11 + n12; n_1 = n11 + n21; n = n11 + n12 + n21 + n22
    upper = min(n1_, n_1)
    lower = max(n1_ + n_1 - n, 0)
    lower == upper && return 1.0
    s = _HyperState(0, 0, 0, 0, 0.0, false)
    prob = _hyper0!(s, n11, n1_, n_1, n)
    sleft = 0.0
    p = _hyper!(s, lower)
    i = lower + 1
    while p <= 0.99999999 * prob
        sleft += p
        p = _hyper!(s, i)
        i += 1
    end
    p <= 1.00000001 * prob && (sleft += p)
    sright = 0.0
    p = _hyper!(s, upper)
    j = upper - 1
    while p <= 0.99999999 * prob
        sright += p
        p = _hyper!(s, j)
        j -= 1
    end
    p <= 1.00000001 * prob && (sright += p)
    return min(sleft + sright, 1.0)
end
