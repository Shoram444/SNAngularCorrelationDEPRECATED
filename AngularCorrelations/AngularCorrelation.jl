using StatsPlots, UnROOT, StatsBase, MPSensitivity, CategoricalArrays, Revise
using FHist, LaTeXStrings, MPThemes, DataFrames, DataFramesMeta, Distributions

include("MiscFuncs.jl")
using .MiscFuncs

plotlyjs()
my_vibrant(;
    size           = (800, 600),
    legend         = :outertopleft,
    guidefontsize  = 8,
    tickfontsize   = 8,
    titlefontsize  = 8,
    legendfontsize = 6,
    left_margin    = 4Plots.mm,
    right_margin   = 4Plots.mm,
    top_margin     = 4Plots.mm,
    bottom_margin  = 4Plots.mm,
    dpi            = 200,
);




f = ROOTFile(
    "/home/shoram/Work/PhD_Thesis/Job15/AngularCorrelations/AngularCorrelationAllEnergies96MilEvents.root",
);
tree = DataFrame(LazyTree(f, "tree", keys(f["tree"])));
@transform! tree :ESum = :reconstructedEnergy2 + :reconstructedEnergy1;

dEmitted = 1 # dθdif in degrees
nBins    = Int(180 / dEmitted)
minAngle = 0
maxAngle = 180
binWidth = maxAngle / nBins

minEnergy = 500
maxEnergy = 3500
dEnergy   = 500

colors = [palette(:seaborn_bright)[i] for i in 1:length(palette(:seaborn_bright))]

###############################################################
################# COMPUTE individual stats

df_stats = DataFrame(;
    minThetaEmitted = Real[],
    maxThetaEmitted = Real[],
    minE = Real[],
    maxE = Real[],
    mean = Real[],
    mode = Real[],
    median = Real[],
    variance = Real[],
    escEfficiency = Real[],
)

for minThetaEmitted in minAngle:dEmitted:(maxAngle - dEmitted)  # for loop over emitted slices, the size of slize is determined by dEmitted 
    maxThetaEmitted = minThetaEmitted + dEmitted

    for minE in minEnergy:dEnergy:(maxEnergy - dEnergy)     # for loop over energy slices of the sum of energies 
        maxE = minE + dEnergy

        push!(
            df_stats,
            get_slice_stats(minThetaEmitted, maxThetaEmitted, minE, maxE, tree, binWidth),
        )
    end
end


###############################################################
################# DRAW gs

dθ = 1
sign = "p"
maxSteps = Int(180 / dθ)
labels = []

fh2ds = Array{Hist2D, 1}(undef, maxSteps)
gs = []
ks = []

gs_cdf = []

gs_abs = []
ks_abs = []

for (i, n) in enumerate(1:dθ:180)
    @show cutEdges1 = get_cut_edges(n - 1, 1, dθ, sign)

    sdf = @chain tree begin  # filter out the dataframe
        @subset((cutEdges1[1] .<= :thetaEscaped .<= cutEdges1[2]))
        @select(:thetaEscaped, :thetaEmitted)
    end

    push!(labels, string(
        "θ ∈ (",
        cutEdges1[1],
        ", ",
        cutEdges1[2],
        " )",
    ))

    fh2ds[i] = Hist2D(
        (sdf[!, :thetaEmitted], sdf[!, :thetaEscaped]),
        (minAngle:dEmitted:maxAngle, minAngle:dEmitted:maxAngle),
    )

    push!(gs, get_diagonal_sums(fh2ds[i]))
    push!(ks, get_k_factors(fh2ds[i]))
    push!(gs_cdf, get_diagonal_sums_cdf(fh2ds[i]))
    push!(gs_abs, get_diagonal_sums_abs(fh2ds[i]))
    push!(ks_abs, get_k_factors_abs(fh2ds[i]))
end


yMax = Int(180 / dθ - 1)
x = -179:1:179
y = 0:yMax

z = zeros(length(y), 359)
for c in eachindex(x)
    for r in eachindex(y)
        z[r, c] = MiscFuncs.get_gs(y[r], x[c], gs)
    end
end

contourStep = 0.1
z1 = map(x -> abs(2 * (0.5 - x)), z) # transform cdf to have maximum at 0.5 and move symetrically to the sides

c = plot(
    x,
    y .* dθ,
    z;
    ylims = (0, 180),
    yticks = 0:15:180,
    xlabel = "k",
    ylabel = "ϕ",
    legend = :none,
    title = string("dϕ= ", dθ, "° "),
    dpi = 150,
    linetype = :surface,
    levels = 0:contourStep:1,
    contour_labels = true,
)
vline!([0], label = "", c = :black, l2 = 4)
savefig(
    string(
        "/home/shoram/Work/PhD_Thesis/Job15/AngularCorrelations/LargeStats/gs/3d_contour_surface",
        dθ,
        "deg.svg",
    ),
)


halfPoints = map( x -> abs.(0.5 .- x), gs_cdf ) # move each point by 0.5 and take abs. 
                                                # This leads to array between 0 and 0.5. 
                                                # Where number closest to 0 represents the 50% statistics.
halfPoints = -179 .+ argmin.(halfPoints[:])     # To find the x position of the 50% stats, shift by 179


modTree = @chain tree begin
    @select(:thetaEmitted, :thetaEscaped, :ESum)
    @rtransform :thetaEscapedModified = shift_angle( :thetaEscaped, halfPoints )
    @subset( 0 .<= :thetaEscapedModified .<= 180)
end