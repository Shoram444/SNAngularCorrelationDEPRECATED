using StatsPlots, UnROOT, StatsBase, MPSensitivity, CategoricalArrays, Revise
using FHist, LaTeXStrings, MPThemes, DataFrames, DataFramesMeta, Distributions

include("MiscFuncs.jl")
using .MiscFuncs

gr()
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

xPts = minAngle:dEmitted:(maxAngle - dEmitted) # x-points for plotting


labelsEnergy = [
    string("E ∈ (", minE, ", ", minE + dEnergy, " )") for
    minE in minEnergy:dEnergy:(maxEnergy - dEnergy)
]

labelsAngles = [
    string("θ ∈ (", minA, ", ", minA + dEmitted, " )") for
    minA in minAngle:dEmitted:(maxAngle - dEmitted)
]

colors = [palette(:seaborn_bright)[i] for i in 1:length(palette(:seaborn_bright))]
#colors = [:red, :orange, :yellow, :green, :cyan, :blue, :purple]

h2d = histogram2d(
    tree.thetaEmitted,
    tree.thetaEscaped;
    nbins        = (nBins, nBins),
    xlabel       = "θemitted",
    ylabel       = "θescaped",
    legend       = :topright,
    title        = string("θesc vs θemit, ", nrow(tree), " entries"),
    lims         = (0, 180),
    aspect_ratio = 1,
)
savefig(h2d, joinpath("AngularCorrelations/LargeStats", string("h2d.png")))

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

df_stats


###############################################################
################# DRAW individual stats

# anim = @animate for minThetaEmitted in minAngle:dEmitted:(maxAngle - dEmitted)
for minThetaEmitted in minAngle:dEmitted:(maxAngle - dEmitted)  # for loop over emitted slices, the size of slize is determined by dEmitted 
    maxThetaEmitted = minThetaEmitted + dEmitted

    for minE in minEnergy:dEnergy:(maxEnergy - dEnergy)     # for loop over energy slices of the sum of energies 
        maxE = minE + dEnergy

        theta_esc = get_theta_esc_slice(minThetaEmitted, maxThetaEmitted, minE, maxE, tree)      # get the vertical slice in θemit

        hist = Hist1D(theta_esc, (minAngle:binWidth:maxAngle))   # get 1D histogram from the slice. This is used to obtain counts

        sh = stephist(
            theta_esc;
            nbins  = nBins,
            c      = :blue,
            xlabel = "θEsc",
            ylabel = string("counts/", binWidth, "°"),
            label  = "SimData",
            title  = string("θemitted ∈ (", minThetaEmitted, ", ", maxThetaEmitted, ")", " total events = ", length(theta_esc), " E ∈ (", minE, ", ", maxE, ")"),
            ylims  = (0, 1.1 * maximum(bincounts(hist))),
            xlims  = (minAngle, maxAngle),
        )

        stats = get_slice_stats(
            minThetaEmitted,
            maxThetaEmitted,
            minE,
            maxE,
            theta_esc;
            _binWidth = binWidth,
        )
        meanData = stats[5]
        modeData = stats[6]
        medianData = stats[7]
        varianceData = stats[8]

        vline!([medianData]; c = :black, label = "median", lw = 2)
        vline!([modeData]; c = :green, label = "mode", lw = 2)
        vline!([meanData]; c = :red, label = "mean", lw = 2)

        annotate!([(
            -45,
            0.6 * maximum(bincounts(hist)),
            (
                string(
                    "mean = ",
                    meanData,
                    "\n median = ",
                    medianData,
                    "\n mode = ",
                    modeData,
                    "\n",
                    L"\sqrt{\frac{1}{n-1}\sum^n{(x_i - \bar{x})^2}} = ",
                    varianceData,
                ),
                6,
            ),
        )])

        outFileName = string(
            "hist-",
            minThetaEmitted,
            "_nBins-",
            nBins,
            "_E-",
            minE,
            "_",
            maxE,
            ".png",
        )
        outDir = string(
            "/home/shoram/Work/PhD_Thesis/Job15/AngularCorrelations/Figs_individual_Escapes/",
            "dTheta_",
            dEmitted,
            "-nBins_",
            nBins,
            "_E-",
            minE,
            "_",
            maxE,
        )

        if (isdir(outDir))
            savefig(sh, joinpath(outDir, outFileName))
        else
            mkdir(joinpath(pwd(), outDir))
            savefig(sh, joinpath(outDir, outFileName))
        end
    end
end

###############################################################
################# DRAW simple stats

plotMeans = [
    df_stats[(df_stats.minE .== minE), :mean] for
    minE in minEnergy:dEnergy:(maxEnergy - dEnergy)
] # creates a matrix of n rows, where each row represents the energy slice
plotModes = [
    df_stats[(df_stats.minE .== minE), :mode] for
    minE in minEnergy:dEnergy:(maxEnergy - dEnergy)
]
plotMedians = [
    df_stats[(df_stats.minE .== minE), :median] for
    minE in minEnergy:dEnergy:(maxEnergy - dEnergy)
]
plotVars = [
    df_stats[(df_stats.minE .== minE), :variance] for
    minE in minEnergy:dEnergy:(maxEnergy - dEnergy)
]
plotEffs = [
    df_stats[(df_stats.minE .== minE), :escEfficiency] for
    minE in minEnergy:dEnergy:(maxEnergy - dEnergy)
]

residualMeans = [(plotMeans[i] .- y(xPts)) ./ plotVars[i] for i in eachindex(plotMeans)] # residuals in number of sigmas
residualModes = [(plotModes[i] .- y(xPts)) ./ plotVars[i] for i in eachindex(plotModes)]
residualMedians =
    [(plotMedians[i] .- y(xPts)) ./ plotVars[i] for i in eachindex(plotMedians)]

l = @layout [
    a{0.7h}
    b
]

mean1 = scatter(
    minAngle:dEmitted:(maxAngle - dEmitted),
    plotMeans;
    label  = reshape(labelsEnergy, 1, length(labelsEnergy)),  # labels must be in matrix form -> hence reshape
    lims   = (minAngle, maxAngle),
    ylabel = L"\theta_{escaped}",
    title  = L"\mathrm{means~vs~energy}",
    legend = :topleft,
    c      = reshape(colors, 1, length(colors)),
    ms     = 2,
)
plot!(y; lw = 2, c = :black, label = L"f(x) = x");

mean2 = plot(
    minAngle:dEmitted:(maxAngle - dEmitted),
    residualMeans;
    xlims = (minAngle, maxAngle),
    ylabel = L"(y_i - f(x_i))/\sigma",
    label = "",
    legend = :topleft,
    xlabel = L"\theta_{emitted}",
    c = reshape(colors, 1, length(colors)),
    seriestype = :steppre,
    lw = 2,
)
hline!([0]; c = :black, lw = 2, label = "")
hline!(
    [1, -1];
    c = :black,
    lw = 1,
    ls = :dash,
    fill = 0,
    fillcolor = :black,
    fillalpha = 0.1,
    label = "",
);

pMean = plot(mean1, mean2; layout = l)

savefig(pMean, "AngularCorrelations/pMean.png")

mode1 = scatter(
    minAngle:dEmitted:(maxAngle - dEmitted),
    plotModes;
    label = reshape(labels, 1, length(plotModes)),                                  # labels must be in matrix form -> hence reshape
    lims = (minAngle, maxAngle),
    ylabel = L"\theta_{escaped}",
    title = L"\mathrm{modes~vs~energy}",
    legend = :topleft,
    c = reshape(colors, 1, length(colors)),
    ms = 2,
);
plot!(y; lw = 2, c = :black, label = L"f(x) = x");

mode2 = plot(
    minAngle:dEmitted:(maxAngle - dEmitted),
    residualModes;
    xlims = (minAngle, maxAngle),
    ylabel = L"(y_i - f(x_i))/\sigma",
    label = "",
    legend = :topleft,
    xlabel = L"\theta_{emitted}",
    c = reshape(colors, 1, length(colors)),
    seriestype = :steppre,
);
hline!([0]; c = :black, lw = 2, label = "");
hline!(
    [1, -1];
    c = :black,
    lw = 1,
    ls = :dash,
    fill = 0,
    fillcolor = :black,
    fillalpha = 0.1,
    label = "",
);

pMode = plot(mode1, mode2; layout = l);

savefig(pMode, "AngularCorrelations/pMode.png")

median1 = scatter(
    minAngle:dEmitted:(maxAngle - dEmitted),
    plotMedians;
    label = reshape(labels, 1, length(plotMedians)),                                  # labels must be in matrix form -> hence reshape
    lims = (minAngle, maxAngle),
    ylabel = L"\theta_{escaped}",
    title = L"\mathrm{medians~vs~energy}",
    legend = :topleft,
    c = reshape(colors, 1, length(colors)),
    ms = 2,
);
plot!(y; lw = 2, c = :black, label = L"f(x) = x");

median2 = plot(
    minAngle:dEmitted:(maxAngle - dEmitted),
    residualMedians;
    xlims = (minAngle, maxAngle),
    ylabel = L"(y_i - f(x_i))/\sigma",
    label = "",
    legend = :topleft,
    xlabel = L"\theta_{emitted}",
    c = reshape(colors, 1, length(colors)),
    seriestype = :steppre,
)
hline!([0]; c = :black, lw = 2, label = "");
hline!(
    [1, -1];
    c = :black,
    lw = 1,
    ls = :dash,
    fill = 0,
    fillcolor = :black,
    fillalpha = 0.1,
    label = "",
);

pMedian = plot(median1, median2; layout = l);

savefig(pMedian, "AngularCorrelations/pMedian.png")

###############################################################
################# DRAW energy histograms

fh1Emitted = Hist1D(tree.thetaEmitted, (minAngle:maxAngle))
h1Emitted = stephist(
    tree.thetaEmitted;
    nbins = nBins,
    xlims = (minAngle, maxAngle),
    ylims = (0, 1.1 * maximum(bincounts(fh1Emitted))),
    xlabel = "θemitted",
    ylabel = string("counts/", binWidth, "°"),
    title = string("θemitted ; events = ", nrow(tree)),
    legend = :false,
    c = :blue,
    lw = 2,
)
savefig(h1Emitted, "AngularCorrelations/h1Emitted.png")

fh1Escaped = Hist1D(tree.thetaEscaped, (minAngle:maxAngle))
h1Escaped = stephist(
    tree.thetaEscaped;
    nbins = nBins,
    xlims = (minAngle, maxAngle),
    ylims = (0, 1.1 * maximum(bincounts(fh1Escaped))),
    xlabel = "θEscaped",
    ylabel = string("counts/", binWidth, "°"),
    title = string("θescaped ; events = ", nrow(tree)),
    legend = :false,
    c = :red,
    lw = 2,
)
savefig(h1Escaped, "AngularCorrelations/h1Escaped.png")

h2dEnergy = histogram2d(
    tree.reconstructedEnergy1,
    tree.reconstructedEnergy2;
    nbins = (350, 350), # 10kev/bin
    lims = (0, 3500),
    xlabel = "E1 [keV]",
    ylabel = "E2 [keV]",
    title = "CD: single-electron spectrum",
    legend = :topright,
    aspect_ratio = 1,
)

savefig(h2dEnergy, "AngularCorrelations/h2dEnergy.png")

###############################################################
################# DRAW efficiencies of cuts

dAngle = 5 # dθdif in degrees
effs = []
maxeff = 0

for (i, minE) in enumerate(minEnergy:dEnergy:(maxEnergy - dEnergy))     # for loop over energy slices of the sum of energies 
    maxE = minE + dEnergy

    eff = Float64[]
    for minThetaEmitted in minAngle:dAngle:(maxAngle - dAngle)  # for loop over emitted slices, the size of slize is determined by dAngle 
        maxThetaEmitted = minThetaEmitted + dAngle

        ε = get_reco_efficiency(minThetaEmitted, maxThetaEmitted, minE, maxE, tree)
        if (maxeff < ε && i != 1)
            maxeff = ε
        end
        push!(eff, ε)
    end
    push!(effs, eff)

end

pThetaEff = plot(
    minAngle:dAngle:(maxAngle - dAngle),
    effs;
    seriestype = :steppre,
    xlabel     = "θemitted",
    ylabel     = "ε",
    title      = string("efficiency of θescaped being within ", dAngle, "° bins of θemitted"),
    label      = reshape(labels, 1, length(labels)),
    legend     = :topright,
    xlims      = (minAngle, maxAngle),
    ylims      = (0, 1.6 * maxeff),
)
savefig(pThetaEff, "AngularCorrelations/pThetaEff.png")

###############################################################
################# DRAW h2ds for various energy ranges


for (i, minE) in enumerate(minEnergy:dEnergy:(maxEnergy - dEnergy))     # for loop over energy slices of the sum of energies 
    maxE = minE + dEnergy
    h2d  = histogram2d(tree[(minE .< tree.ESum .< maxE), 5], tree[(minE .< tree.ESum .< maxE), 10]; nbins = (nBins, nBins), xlabel = "θemitted", ylabel = "θescaped", legend = :topright, title = string("θesc vs θemit, ", " entries", length(tree[(minE .< tree.ESum .< maxE), 5]), "; E ∈", minE, " - ", maxE, " keV"), lims = (minAngle, maxAngle), aspect_ratio = 1)

    outFileName = string("h2d-", minE, "-", maxE, ".png")
    outDir      = string("./AngularCorrelations/LargeStats/h2d/")

    plot!(y; lw = 2, c = :black, label = "")

    if (isdir(outDir))
        savefig(h2d, joinpath(outDir, outFileName))
    else
        mkdir(outDir)
        savefig(h2d, joinpath(outDir, outFileName))
    end
end


###############################################################
################# DRAW directions each-side, same-side of foil!

zs = zeros(nrow(tree))

h2dAngles = Hist2D((tree.thetaEmitted, tree.thetaEscaped), (0:1:180, 0:1:180))
histogram2d(tree.thetaEmitted, tree.thetaEscaped; lims = (0, 180))

h2dDirsAbsolute = Hist2D(
    (zs, zs), # initialize empty histogram
    (0:1:180, 0:1:180),
)
bincounts(h2dDirsAbsolute)[1, 1] = 0

for row in eachrow(tree)
    if (row.momentumEscaped1x * row.momentumEscaped2x > 0)
        push!(h2dDirsAbsolute, row.thetaEmitted, row.thetaEscaped, 1)
    elseif (row.momentumEscaped1x * row.momentumEscaped2x < 0)
        push!(h2dDirsAbsolute, row.thetaEmitted, row.thetaEscaped, 0)
    else
        continue
    end
end

h2dDirsAbsolute

h2dDirectionsAbsolute = plot(
    h2dDirsAbsolute;
    xlabel = "θemitted",
    ylabel = "θescaped",
    title = "escape on each side of foil: 0, on same side: +1  ",
    lims = (0, 180),
    legend = :topright,
    aspect_ratio = 1,
)
plot!(y; lw = 2, c = :black, label = "")
savefig(h2dDirectionsAbsolute, "AngularCorrelations/LargeStats/h2dDirectionsAbsolute.png")


h2dDirsRelative = Hist2D(
    (zs, zs), # initialize empty histogram
    (0:1:180, 0:1:180),
)
bincounts(h2dDirsRelative)[1, 1] = 0.0

for row in eachrow(tree)
    if (row.momentumEscaped1x * row.momentumEscaped2x > 0)
        push!(h2dDirsRelative, row.thetaEmitted, row.thetaEscaped, 1)
    elseif (row.momentumEscaped1x * row.momentumEscaped2x < 0)
        push!(h2dDirsRelative, row.thetaEmitted, row.thetaEscaped, 0)
    else
        continue
    end
end

h2dDirectionsRelative = plot(
    h2dDirsRelative;
    xlabel = "θemitted",
    ylabel = "θescaped",
    title = "escape on each side of foil: 0, on same side: 1  ",
    lims = (0, 180),
    legend = :topright,
    aspect_ratio = 1,
)
plot!(y; lw = 2, c = :black, label = "")
savefig(h2dDirectionsRelative, "AngularCorrelations/LargeStats/h2dDirectionsRelative.png")

mat = convert.(Float64, bincounts(h2dDirsRelative))

for x in 1:1:180, y in 1:1:180
    rel = mat[Int(x), Int(y)]
    n = bincounts(h2dAngles)[Int(x), Int(y)]
    mat[Int(x), Int(y)] = rel / n
end


heatmap(
    mat';
    clim         = (0, 1),
    lims         = (0, 180),
    xlabel       = "θemitted",
    ylabel       = "θescaped",
    title        = "escape on each side of foil: 0, on same side: +1  ",
    legend       = :topright,
    aspect_ratio = 1,
)
plot!(y; lw = 2, c = :black, label = "")
savefig("Allangles.png")

bincounts(h2dAngles)

subTree_30deg = @subset(
    tree,
    abs.(
        :momentumEscaped1x ./
        sqrt.(
            :momentumEscaped1x .^ 2 .+ :momentumEscaped1y .^ 2 .+ :momentumEscaped1z .^ 2
        )
    ) .> cos(deg2rad(30)),
    abs.(
        :momentumEscaped2x ./
        sqrt.(
            :momentumEscaped2x .^ 2 .+ :momentumEscaped2y .^ 2 .+ :momentumEscaped2z .^ 2
        )
    ) .> cos(deg2rad(30)),
)


h2dDirsRelative30 = Hist2D(
    (zs, zs), # initialize empty histogram
    (0:1:180, 0:1:180),
)
bincounts(h2dDirsRelative30)[1, 1] = 0.0

for row in eachrow(subTree_30deg)
    if (row.momentumEscaped1x * row.momentumEscaped2x > 0)
        push!(h2dDirsRelative30, row.thetaEmitted, row.thetaEscaped, 1)
    elseif (row.momentumEscaped1x * row.momentumEscaped2x < 0)
        push!(h2dDirsRelative30, row.thetaEmitted, row.thetaEscaped, 0)
    else
        continue
    end
end

mat = convert.(Float64, bincounts(h2dDirsRelative30))
h2dAngles =
    Hist2D((subTree_30deg.thetaEmitted, subTree_30deg.thetaEscaped), (0:1:180, 0:1:180))
for x in 1:1:180, y in 1:1:180
    rel = mat[Int(x), Int(y)]
    n = bincounts(h2dAngles)[Int(x), Int(y)]
    mat[Int(x), Int(y)] = rel / n
end


heatmap(
    mat';
    clim         = (0, 1),
    lims         = (0, 180),
    xlabel       = "θemitted",
    ylabel       = "θescaped",
    title        = "escape on each side of foil: 0, on same side: +1  ",
    legend       = :topright,
    aspect_ratio = 1,
)
plot!(y; lw = 2, c = :black, label = "")
savefig("FilterecCos30.png")

###############################################################
################# DRAW gs

midPoint1 = 0   # starting point for the angle cut ( θ = mid ± k*dθ)
# midPoint2 = 150
dθ = 6
sign = "p"
maxSteps = 180 #Int(midPoint1 / dθ) 
labels = []

fh2ds = Array{Hist2D, 1}(undef, maxSteps)
gs = []
ks = []

gs_abs = []
ks_abs = []

@time for n in 1:dθ:maxSteps
    @show cutEdges1 = get_cut_edges(n - 1, 1, dθ, sign)
    # cutEdges2 = get_cut_edges(midPoint2, n, dθ, sign)

    sdf = @chain tree begin  # filter out the dataframe
        @subset((cutEdges1[1] .<= :thetaEscaped .<= cutEdges1[2]))#.| 
        # (cutEdges2[1] .<= :thetaEscaped .<= cutEdges2[2]),)
        @select(:thetaEscaped, :thetaEmitted)
    end

    push!(labels, string(
        "θ ∈ (",
        cutEdges1[1],
        ", ",
        cutEdges1[2],
        " )",
        # " + (",
        # cutEdges2[1],
        # ", ",
        # cutEdges2[2],
        # " )",
    ))

    fh2ds[n] = Hist2D(
        (sdf[!, :thetaEmitted], sdf[!, :thetaEscaped]),
        (minAngle:dEmitted:maxAngle, minAngle:dEmitted:maxAngle),
    )

    push!(gs, get_diagonal_sums(fh2ds[n]))
    push!(ks, get_k_factors(fh2ds[n]))
    push!(gs_abs, get_diagonal_sums_abs(fh2ds[n]))
    push!(ks_abs, get_k_factors_abs(fh2ds[n]))
end

# plot(
#     ks,
#     gs;
#     label = reshape(labels, 1, length(labels)),
#     palette = :matter,
#     # c = reshape(colors, 1, length(colors)),
#     xlims = (-179, 179),
#     xlabel = "k-factor",
#     ylabel = "g(k)",
#     title = "Sum over diagonals of the f(θ,ϕ) function ",
#     legend = :topright,
#     size = (3600, 2600),
#     dpi = 100,
#     lw = 10,
#     alpha = 0.4,
# )

# t = string(midPoint1)#, "_and_", midPoint2)

# savefig(
#     string(
#         "/home/shoram/Work/PhD_Thesis/Job15/AngularCorrelations/LargeStats/gs/gs_angle",
#         t,
#         "_p1.png",
#     ),
# )



# p1 = plot()
# for i in 1:length(ks_abs)
#     plot!(
#         ks_abs[i],
#         gs_abs[i] / maximum(gs_abs[i]);
#         label = labels[i],
#         c = palette(:matter)[i],
#         # c = colors[i],
#         # seriestype = :steppre,
#         xlims = (0, 179),
#         xlabel = "|k-factor|",
#         ylabel = "g(k)",
#         title = "Sum over diagonals of the f(θ,ϕ) function ",
#         legend = :topright,
#         size = (3600, 2400),
#         dpi = 400,
#         lw = 3,
#         alpha = 0.4,
#     )
# end
# p2 = scatter(
#     [0, 0],
#     [0, 1],
#     zcolor = [0, 180],
#     clims = (0, 180),
#     xlims = (1, 1.1),
#     # xshowaxis = false,
#     # yshowaxis = false,
#     framestyle = :none,
#     label = "",
#     c = :matter,
#     grid = false,
# )
# l = @layout [a b{0.01w}]
# p = plot(p1, p2, layout = l, size = (3600, 2400), dpi = 100)
# savefig(
#     p,
#     string(
#         "/home/shoram/Work/PhD_Thesis/Job15/AngularCorrelations/LargeStats/gs/gs_ABS_angle_",
#         t,
#         "_p1.png",
#     ),
# )



yMax = Int(180 / dθ - 1)
x = -179:1:179
y = 0:yMax

z1 = zeros(length(y), 359)
for c in eachindex(x)
    for r in eachindex(y)
        z1[r, c] = get_gs_rel(y[r], x[c], gs)
    end
end

contour(
    x,
    y .* dθ,
    z1;
    ylims = (0,180),
    yticks = 0:15:180,
    xlabel = "k",
    ylabel = "ϕ",
    legend = :right,
    colorbar_title = "g(k, ϕ)",
    dpi = 150,
)
vline!([0], label = "", c = :black, l2 = 4)
savefig(
    string(
        "/home/shoram/Work/PhD_Thesis/Job15/AngularCorrelations/LargeStats/gs/contour_",
        dθ,
        "deg.png",
    ),
)

