module MiscFuncs
using DataFrames, DataFramesMeta, FHist, StatsBase

export y,
    get_theta_esc_slice,
    Gaus,
    get_chi2,
    get_reco_efficiency,
    get_diagonal_sums,
    get_k_factors,
    get_diagonal_sums_abs,
    get_k_factors_abs,
    get_slice_stats,
    get_cut_edges,
    get_gs


y(x) = x # generic linear line y = x 

function get_theta_esc_slice(_θemiMin::Real, _θemiMax::Real, _tree::DataFrame)
    return θesc = _tree[
        (_tree.thetaEmitted .> _θemiMin) .& (_tree.thetaEmitted .< _θemiMax),
        :thetaEscaped,
    ]
end
function get_theta_esc_slice(
    _θemiMin::Real,
    _θemiMax::Real,
    _EMin::Real,
    _EMax::Real,
    _tree::DataFrame,
)
    return θesc = _tree[
        (_tree.thetaEmitted .> _θemiMin) .& (_tree.thetaEmitted .< _θemiMax) .& (_tree.ESum .> _EMin) .& (_tree.ESum .< _EMax),
        :thetaEscaped,
    ] # return slice where all conditions are fulfilled
end

function Gaus(fit, x, N)
    return N / (fit.σ * sqrt(2 * π)) * exp(-0.5 * ((x - fit.μ) / fit.σ)^2)
end

function Gaus(μ, σ, x, N)
    return N / (σ * sqrt(2 * π)) * exp(-0.5 * ((x - μ) / σ)^2)
end

function get_chi2(x::Vector{<:Real}, y::Vector{<:Real}, f::Function)
    xi = 0
    ndf = 0
    for idx in eachindex(y)
        if (y[idx] <= 0)
            continue
        else
            xi += (y[idx] - f(x[idx]))^2 / y[idx]
            ndf += 1
        end
    end
    return xi, ndf - 3 # there are 3 fit parameters (μ. σ, N) 
end

function get_chi2(y::Vector{<:Real}, fx::Vector{<:Real})
    xi = 0
    ndf = 0
    for idx in eachindex(y)
        if (y[idx] <= 0)
            continue
        else
            xi += (y[idx] - fx[idx])^2 / y[idx]
            ndf += 1
        end
    end
    return xi, ndf - 3 # there are 3 fit parameters (μ. σ, N) 
end

"""
get_pmSigma_mode(_hist; _sigmaHalf=0.3414)

## Description of get_pmSigma_mode

Returns how many bins to the left and to the right from mode it is to cover ±σ (σ = 0.3414, or as stated in argument).
The returned values are [l,r] vector where l represents mode-sigmaBin, and r represents sigmaBin-mode.

Input arguments are:

  - Hist1D objed from FHist.jl
  - Sigma portion (by default 0.3414)
"""
function get_pmSigma_mode(_hist; _sigmaHalf = 0.3414)
    totalEvents = sum(bincounts(_hist))
    mode_bin = 0 # in which bin the mode is placed

    if (length(argmax(bincounts(_hist))) == 1 || sum(bincounts(_hist)) == 0) # makes sure there is a mode to begin with, if not return 0
        mode_bin = argmax(bincounts(_hist))
    else
        return [0, 0]
    end

    tempSumLeft = (bincounts(_hist)[mode_bin] / 2) / totalEvents   # to begin with the sum is initialized to 1/2 of mode's bin
    tempSumRight = (bincounts(_hist)[mode_bin] / 2) / totalEvents

    binLeftCounter = mode_bin - 1
    binRightCounter = mode_bin + 1

    while (tempSumLeft <= _sigmaHalf && binLeftCounter != 0) #68.27% / 2 half of one sigma
        tempSumLeft += bincounts(_hist)[binLeftCounter] / total # add the content of the bin to the left
        binLeftCounter -= 1
    end

    while (tempSumRight <= _sigmaHalf && binRightCounter != length(bincounts(_hist))) #68.27% / 2 half of one sigma
        tempSumRight += bincounts(_hist)[binRightCounter] / total # add the content of the bin to the right
        binRightCounter += 1
    end

    return [mode_bin - binLeftCounter, binRightCounter - mode_bin]
end

"""
get_reco_efficiency(_θemiMin::Real, _θemiMax::Real,_slice::Vector{Real})

## Description of get_reco_efficiency

Returns the efficiency of the reconstructed slice as:

`ε = (n of events escaped within (θ_{min}, θ_{max})/(n of all events in (θ_{min}, θ_{max})`

Input arguments are:

  - Hist1D objed from FHist.jl
  - Sigma portion (by default 0.3414)
"""
function get_reco_efficiency(_θemiMin::Real, _θemiMax::Real, _slice::Vector{<:Real})
    n_passed = filter(θ -> (_θemiMin <= θ <= _θemiMax), _slice) |> length
    return length(_slice) > 0 ? n_passed / length(_slice) : 0
end

function get_reco_efficiency(
    _θemiMin::Real,
    _θemiMax::Real,
    _EMin::Real,
    _EMax::Real,
    _tree::DataFrame,
)
    return get_reco_efficiency(
        _θemiMin,
        _θemiMax,
        get_theta_esc_slice(_θemiMin, _θemiMax, _EMin, _EMax, _tree),
    )
end

"""
get_diagonal_sums(h2d)

## Description of get_diagonal_sums

Returns a vector of the sums of the diagonals of the 2d histogram.
(Diagonals from top-left to bottom-right)

Input arguments are:

  - Hist2D object from FHist.jl
"""
function get_diagonal_sums(h2d)
    mat = bincounts(h2d)      # get the matrix of bin counts of the 2d histogram
    m = size(mat)[1]        # get the number of rows/cols of mxm matrix
    ndiagonals = 2 * (m) - 1       # there are 2m-1 diagonals in a mxm matrix
    sums = zeros(ndiagonals)   # initialize the array of sums of the diagonal elements

    for r in 1:m           # loop over rows r
        for c in 1:m       # loop over columns c
            sums[m - (r - c)] += mat[r, c]
        end
    end
    return sums
end

function get_diagonal_sums_abs(h2d)
    mat = bincounts(h2d)       # get the matrix of bin counts of the 2d histogram
    m = size(mat)[1]         # get the number of rows/cols of mxm matrix
    sums = zeros(m + 1)           # initialize the array of sums of the diagonal elements
    # there are m+1 diagonals in a mxm matrix

    for i in 1:m                # loop over rows i
        for j in 1:m            # loop over columns j
            sums[Int(abs(i - j) + 1)] += mat[i, j]
        end
    end
    return sums
end

"""
get_k_factors(h2d)

## Description of get_k_factors

Returns a vector of the k-factors, representing by how much there is over/under-estimate.
(Diagonals from top-left to bottom-right)

Input arguments are:

  - Hist2D object from FHist.jl
"""
function get_k_factors(h2d)
    mat = bincounts(h2d)      # get the matrix of bin counts of the 2d histogram
    m = size(mat)[1]        # get the number of rows/cols of mxm matrix
    ks = zeros(2 * (m) - 1)      # initialize array of k-factors

    for d in 1:(2m - 1)            # loop over the diagnoal index in sums
        k = m - d
        ks[d] = k
    end
    return ks
end

function get_k_factors_abs(h2d)
    mat = bincounts(h2d)      # get the matrix of bin counts of the 2d histogram
    m = size(mat)[1]        # get the number of rows/cols of mxm matrix

    return [k for k in 0:m]
end

function get_slice_stats(
    _θemiMin::Real,
    _θemiMax::Real,
    _EMin::Real,
    _EMax::Real,
    _theta_esc::Vector{Float64},
    _binWidth = 1.0,
)
    meanData = 0
    modeData = 0
    medianData = 0
    varianceData = 0
    escEffData = 0

    hist = Hist1D(_theta_esc, (0:_binWidth:180)) # used to get mode from the binned histogram

    if (length(_theta_esc) > 0)
        meanData = mean(_theta_esc)
        medianData = median(_theta_esc)
        varianceData = sqrt(var(_theta_esc; corrected = :true))
        escEffData = get_reco_efficiency(_θemiMin, _θemiMax, _theta_esc)
        modeData = argmax(bincounts(hist)) * _binWidth
    end

    return [
        _θemiMin,
        _θemiMax,
        _EMin,
        _EMax,
        meanData,
        modeData,
        medianData,
        varianceData,
        escEffData,
    ]
end

function get_slice_stats(
    _θemiMin::Real,
    _θemiMax::Real,
    _EMin::Real,
    _EMax::Real,
    _tree::DataFrame,
    _binWidth = 1.0,
)
    meanData = 0
    modeData = 0
    medianData = 0
    varianceData = 0
    escEffData = 0

    theta_esc = get_theta_esc_slice(_θemiMin, _θemiMax, _EMin, _EMax, _tree)
    hist = Hist1D(theta_esc, (0:_binWidth:180)) # used to get mode from the binned histogram

    if (length(theta_esc) > 0)
        meanData = mean(theta_esc)
        medianData = median(theta_esc)
        varianceData = sqrt(var(theta_esc; corrected = :true))
        escEffData = get_reco_efficiency(_θemiMin, _θemiMax, theta_esc)
        modeData = argmax(bincounts(hist)) * _binWidth
    end

    return [
        _θemiMin,
        _θemiMax,
        _EMin,
        _EMax,
        meanData,
        modeData,
        medianData,
        varianceData,
        escEffData,
    ]
end

function get_cut_edges(midPoint, multiple, dθ, sign)
    l = midPoint - multiple * dθ
    u = midPoint + multiple * dθ
    if (l < 0)
        l = 0
    end
    if (u > 180)
        u = 180
    end

    if (sign == "pm")
        return (l, u)
    elseif (sign == "p")
        return (midPoint, u)
    elseif (sign == "m")
        return (l, midPoint)
    else
        error("Must chose a valid sign from options: pm, m, p")
        return 0
    end
end

function get_gs(θ, k, gs) # return the value of g corresponding to the θ and k
    r = Int(floor(θ))+1
    c = (k+length(gs[1])) - 179
    return gs[r][c]
end

end #MODULE END