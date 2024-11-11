using Plots

gr(label="",size=(450,400), grid=false, colorbar_titlefontsize=15, legendfontsize=15,
    guidefontsize=15, tickfontsize=15,colorbar_tickfontsize=12, linewidth=2,
    yticks=(30:30:60,["$i Mb" for i in 29:30]),
    xticks=(0:30:60, ["$i Mb" for i in 28:1:30]),     
    foreground_color_legend = nothing, background_color_legend=nothing,
    framestyle=:box)

function hic_plot(hic_map, file_name)
    heatmap(hic_map, clims=(0,0.25), color=cgrad(:dense), 
        aspect_ratio=1,colorbartitle="\nContact probability")
    png(file_name)
end

function triplet_plot(triplets, file_name, bait_point; max_c=0.1)
    norm_triplets=triplets[:,:,bait_point]
    norm_triplets=half_half(zeros(size(norm_triplets)), norm_triplets)
    heatmap(norm_triplets, clims=(0, max_c),aspect_ratio=1,
        color=cgrad(:gist_heat, rev=true),
        colorbartitle="\nTriplet probability")
    png(file_name)
end

function average_2D_plot(triplets, file_name)
    heatmap(project_2d(triplets, periodic=true),
        ticks=:auto, colorbartitle="\nTriplet probability",size=(500,400),yrotation=0,
        color=cgrad(:gist_heat, rev=true), clims=(0,0.1), scale=:log10)
    png(file_name)
end

function p_value_plot(p_values, file_name)
    heatmap(p_values, clims=(0,0.15), color=cgrad(:hot, scale=:log10),
        aspect_ratio=1,colorbartitle="\nP-value")
    png(file_name)
end

function bh_p_value_plot(p_values, file_name)
    heatmap(adjust_array(p_values,BenjaminiHochberg()),
        aspect_ratio=1,colorbartitle="\nBenjamini-Hochberg adj. p-value",
        color=cgrad(:tempo,rev=true, scale=:log10), clims=(0,0.3))
    png(file_name)
end

function z_score_plot(zs, bait_point,file_name)
    heatmap(zs[:,:,bait_point], colorbartitle="\nZ-score",
        aspect_ratio=1,color=cgrad(:coolwarm), clims=(-4,4))
    png(file_name)
end

function average_2D_z_score_plot(zs, file_name)
    heatmap(project_2d(zs), scale=:log10,size=(500,400),
        ticks=:auto, colorbartitle="\nMean z-score", yrotation=0,
        color=cgrad(:coolwarm), clims=(-4,4))
    png(file_name)
end

function plot_p_3_s_curves(P_3_s, labels, file_name)
    N=size(P_3_s)[1]
    bin_size=2/N
    xs=(1:N).*bin_size #axis in Mb

    plot(xs,P_3_s, label=labels,
        xlabel="Genomic length largest loop (Mb)", scale=:log10, ticks=:auto,
        ylabel="Mean triplet probability", size=[540,400], yrotation=0,
        legend=:bottomleft)
    Plots.pdf(file_name)
end