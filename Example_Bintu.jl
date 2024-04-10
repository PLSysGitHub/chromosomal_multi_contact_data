"""
Example script for generating plots for predictions of contact triplets
on a linear (section of a) chromosome.

"""

include("ContactTripletPredictions.jl")

using Plots, LaTeXStrings
pythonplot(label="",size=(550,500), grid=false, colorbar_titlefontsize=15, legendfontsize=15,
        guidefontsize=15, tickfontsize=15,colorbar_tickfontsize=12, linewidth=1.5,
        xticks=(0:30:60,28:1:30), yticks=(30:30:60, 29:30), #units Mb
        foreground_color_legend = nothing, background_color_legend=nothing,
        framestyle=:box)

#Change names of files and directories to analyse different data
contact_file="./Contact_files/bintu_IMR90.txt"
P0_file="./Contact_files/P0_linear_n10_b119.27064620249611_N65_R7.021446824210077.txt"
triplet_file="./Triplet_files/bintu_IMR90.h5" #triplet file should contain 3D triplet frequency array as "triplets"
out_dir="./Output/Bintu_example/"
if !isdir(out_dir) mkdir(out_dir) end

#Point around which to build plots
bait_point=40

#Read in the data
P=readdlm(contact_file)
P0=readdlm(P0_file)
triplets=read(h5open(triplet_file,"r"),"triplets")
num_samp=h5read(triplet_file,"num_samples")
f_factors=h5read(triplet_file,"f_factors") #Polovnikov et.al. 2019, correction factor for fractal dimension not 2
f_factors[isnan.(f_factors)].=0
triplets[isnan.(triplets)].=0

#Do all plots and save in out_dir
P_3_s = triplets_1d(triplets, periodic=false)
P_3_s = cat(P_3_s, zeros(length(P_3_s),3), dims=2) #space for predicted P_3(s) curves

heatmap(triplets[:,:,bait_point], clims=(0,0.1),
        color=cgrad(:gist_heat, rev=true),
        colorbartitle="Triplet probability")
png(out_dir*"measured_triplets_bait_$(bait_point)")

#2D average probabilities
heatmap(project_2d(triplets, periodic=false),
        ticks=:auto, colorbartitle="Triplet probability",
        color=cgrad(:gist_heat, rev=true), clims=(0,0.1), scale=:log10)
png(out_dir*"measured_probabilities_2D_average")

for (index,prediction) in enumerate(["ideal", "loop_extr", "pairwise", "quad_hamiltonian"])
        file_prefix=out_dir*prediction
        if prediction== "ideal"
                pred_triplets=ideal(P, periodic=false)
        elseif prediction=="loop_extr"
                pred_triplets=loop_extr(P, periodic=false)
        elseif prediction=="pairwise"
                pred_triplets=pairwise_int(P, P0, periodic=false)
                #save the effective energy map
                heatmap(-log.(P./P0), clims=(-5,5), color=cgrad(:coolwarm))
                png(file_prefix*"_energy_estimate")
        elseif prediction=="quad_hamiltonian"
                pred_triplets=ideal(P, periodic=false)
                pred_triplets.*=f_factors
        else                
                error("Unknown prediction label $prediction")
        end

        #Save the P_3(s) curve for later
        P_3_s[:,index+1]=triplets_1d(pred_triplets, periodic=false)

        #Heatmaps below. Color set for each plot type
        heatmap(pred_triplets[:,:,bait_point], clims=(0,0.1),
                color=cgrad(:gist_heat, rev=true),
                colorbartitle="Triplet probability")
        png(file_prefix*"_triplets_bait_$(bait_point)")
        display(current())

        p_values=p_vals(triplets[:,:,bait_point],pred_triplets[:,:,bait_point],num_samp)

        heatmap(p_values, clims=(0,0.1), color=cgrad(:hot),
        colorbartitle="P-value")
        png(file_prefix*"_p_vals_bait_$(bait_point)")
        display(current())

        heatmap(adjust_array(p_values,BenjaminiHochberg()),
        colorbartitle="Benjamini-Hochberg adj. p-value",
        color=cgrad(:tempo,rev=true), clims=(0,0.3))
        png(file_prefix*"_BH_adj_p_bait_$(bait_point)")
        display(current())

        zs=z_scores(triplets,pred_triplets,num_samp)

        heatmap(zs[:,:,bait_point], colorbartitle="Z-score",
                color=cgrad(:coolwarm), clims=(-4,4))
                png(file_prefix*"_z_scores_bait_$(bait_point)")
        display(current())

        #2D average z-scores
        heatmap(project_2d(zs, periodic=false),
                ticks=:auto, colorbartitle="Mean z-score",
                color=cgrad(:coolwarm), clims=(-4,4))
                png(file_prefix*"_z_scores_2D_average")
        display(current())

        heatmap(project_2d(pred_triplets, periodic=false),
                ticks=:auto, colorbartitle="Triplet probability",
                color=cgrad(:gist_heat, rev=true), clims=(0,0.1), scale=:log10)
        png(file_prefix*"_probabilities_2D_average")
        display(current())
end

#Save all the P_3(s) curves
xs=(4*4:4:4*(size(P_3_s)[1]+3))*10 #axis in kb

plot(xs, P_3_s[:, 1:4], label=["Bintu et.al. data" "Ideal polymer" "Loop-extruder" "Pairwise interaction" "Quadratic Hamiltonian" ],
        xlabel="Genomic length largest loop (kb)", scale=:log10, ticks=:auto,
        ylabel="Mean triplet probability", size=[640,500], legend=:bottomleft)
png(out_dir*"P_3_s_curves")

plot(xs, P_3_s[:,[1,5]], label=["Bintu et.al. data" "Quadratic Hamiltonian" ],
        xlabel="Genomic length largest loop (kb)", scale=:log10, ticks=:auto,
        ylabel="Mean triplet probability", size=[640,500], legend=:bottomleft)
png(out_dir*"quad_hamiltonian_P_3_s_curves")

#Make comparison plots for scalings
P_S=P_s(P, periodic=false)
plot(P_S, scale=:log10, xticks=(1:30:61,0:2), yticks=:auto, 
        ylabel="Bintu et.al. P(s)", xlabel="s (Mb)", size=(500,400))
png(out_dir*"P_s")

#Test P(l_1,l_2) scalings
N=size(P)[1]
counter_example=zeros(Int((N-1)/2),N-1)
for i in 1:N-1
        for j in i:N-1-i
                l_3=i+j
                counter_example[i,j]=P_S[i]*P_S[j]*P_S[l_3]/P0[l_3]
        end
end
counter_example[iszero.(counter_example)].=NaN
heatmap(counter_example, title=L"P(x)P(y)P(s)/P_0(s)", ticks=:auto, 
        color=cgrad(:gist_heat, rev=true), clims=(0,0.16), scale=:log10,
        xlabel=L"x\ \mathrm{(kb)}", ylabel=L"y\ \mathrm{(kb)}",
        xticks=([100/30,1000/30],[L"10^2",L"10^3"]), yticks=([100/30,1000/30],[L"10^2",L"10^3"]))
png(out_dir*"product_averages_2D_pairwise")
pred_triplets=pairwise_int(P, P0, periodic=false)
pred_2d=project_2d(pred_triplets, periodic=false)
pred_2d.*=sum(counter_example[.!isnan.(counter_example)])/sum(pred_2d[.!isnan.(pred_2d)])
heatmap(pred_2d,title="Scaled mean, Pairwise int.", ticks=:auto, 
        color=cgrad(:gist_heat, rev=true), clims=(0,0.16), scale=:log10, 
        xlabel=L"x\ \mathrm{(kb)}", ylabel=L"y\ \mathrm{(kb)}", xticks=([100/30,1000/30],[L"10^2",L"10^3"]), yticks=([100/30,1000/30],[L"10^2",L"10^3"]))
png(out_dir*"scaled_mean_2D_pairwise")
