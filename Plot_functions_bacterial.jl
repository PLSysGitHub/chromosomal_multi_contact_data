using Plots, LaTeXStrings

gr(label="",size=(450,400), grid=false, colorbar_titlefontsize=15, legendfontsize=15,
    guidefontsize=15, tickfontsize=15,colorbar_tickfontsize=11, linewidth=2,
    yticks=(200:200:400,["$i Mb" for i in 2:2:4]),
    xticks=(0:200:400, ["$i Mb" for i in 0:2:4]), yrotation=90,       
    foreground_color_legend = nothing, background_color_legend=nothing,
    framestyle=:box, right_margin = 10Plots.mm)

function hic_plot(hic_map, file_name)
    heatmap(hic_map, clims=(0,0.25), color=cgrad(:dense, scale=:exp), 
        aspect_ratio=1,colorbartitle="\nContact probability")
    png(file_name)
end

function triplet_plot(triplets, file_name, bait_point; max_c=0.005, label="Predicted triplet probability")
    norm_triplets=triplets[:,:,bait_point]
    heatmap(norm_triplets, clims=(0, max_c),aspect_ratio=1, size=(500,400),
        color=cgrad(:gist_heat, rev=true, scale=:exp), colorbar_ticks=(0:0.001:0.005),
        rightmargin=20Plots.mm,colorbartitle="\n\n"*label)
    png(file_name)
end

function average_2D_plot(triplets, file_name)
    heatmap(project_2d(triplets, periodic=true), yrotation=0,
        ticks=:auto, colorbartitle="\nTriplet probability",size=(500,400),
        color=cgrad(:gist_heat, rev=true), clims=(0,0.1), scale=:log10)
    png(file_name)
end

function p_value_plot(p_values, file_name)
    p_value_lower=half_half(p_values, NaN.*ones(size(p_values)))
    heatmap(p_value_lower, clims=(0,0.15), color=cgrad(:hot, scale=:log10),
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
    z_score_lower=half_half(zs[:,:,bait_point], NaN.*ones(size(zs[:,:,bait_point])))
    heatmap(z_score_lower, colorbartitle="\nZ-score",
        aspect_ratio=1,color=cgrad(:coolwarm), clims=(-4,4))
    png(file_name)
end

function average_2D_z_score_plot(zs, file_name)
    heatmap(project_2d(zs), scale=:log10, yrotation=0,
        ticks=:auto, colorbartitle="\nMean z-score", size=(500,400),
        color=cgrad(:coolwarm), clims=(-4,4))
    png(file_name)
end

function plot_p_3_s_curves(P_3_s, labels, file_name)
    xs=(1:length(P_3_s[:,1]))*10 #axis in kb

    plot(xs,P_3_s, label=labels, legend=:bottomleft, linewidth=3, yrotation=0,
        xlabel="Genomic length largest loop (kb)", scale=:log10, ticks=:auto,
        ylabel="Mean triplet probability", size=[540,400])
    Plots.pdf(file_name)
end