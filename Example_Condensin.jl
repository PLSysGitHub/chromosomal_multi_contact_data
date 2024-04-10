"""
Example script for generating plots for predictions of contact triplets
on a bacterial chromosome with condensins.

Script calculates ideal polymer, loop-extruder, and pairwise interaction predictions
Compares each to data, and saves plots in given directory

"""

include("ContactTripletPredictions.jl")
using Plots,HDF5, DelimitedFiles, StatsBase
pythonplot(label="",size=(540,500), grid=false, colorbar_titlefontsize=15, legendfontsize=15,
        guidefontsize=15, tickfontsize=15,colorbar_tickfontsize=12, linewidth=1.5,
        yticks=(200:200:400,["$i kb" for i in 2000:2000:4000]), yrotation=90,
        xticks=(0:200:400, ["$i kb" for i in 0:2000:4000]),#units in kb        
        foreground_color_legend = nothing, background_color_legend=nothing,
        framestyle=:box)

#Change names of files and directories to analyse different data
num_samp=3060
contact_file="./Contact_files/condensin_hic.txt"
P0_file="./Contact_files/P0_condensin.txt"
triplet_file="./Triplet_files/condensin_triplets.h5" #triplet file should contain 3D triplet frequency array as "triplets"
out_dir="./Output/Condensin_example/"
if !isdir(out_dir) mkdir(out_dir) end

#Read in the data
P=readdlm(contact_file)
P0=readdlm(P0_file)
triplets=read(h5open(triplet_file,"r"),"triplets")
f_factors=h5read(triplet_file,"f_factors") #Polovnikov et.al. 2019, correction factor for fractal dimension not 2
f_factors[isnan.(f_factors)].=0
f_factors[isinf.(f_factors)].=0
N=size(P)[1]
#Start pipeline
P_3_s = triplets_1d(triplets)
P_3_s = cat(P_3_s, zeros(length(P_3_s),3), dims=2) #space for predicted P_3(s) curves

#Point around which to build plots
bait_point=300

heatmap(P, clims=(0,0.25), color=cgrad(:dense), colorbartitle="Contact probability", size=(550,400))
png(out_dir*"contact_map")

#Make comparison plots for scalings
P_S=P_s(P, periodic=true)
plot((2:N)*10,P_S, scale=:log10, xticks=:auto, yticks=:auto, ylabel="Condensin simulation P(s)", xlabel="s")
png(out_dir*"P_s")

#Actual triplet probabilities
norm_triplets=triplets[:,:,bait_point]./mean(triplets[:,:,bait_point])
norm_triplets=half_half(norm_triplets, zeros(size(norm_triplets)))
heatmap(norm_triplets, clims=(0,20),
        color=cgrad(:gist_heat, rev=true),
        colorbartitle="Triplet probability relative to mean", size=(550,400))
png(out_dir*"measured_triplets_bait_$(bait_point)")

#2D average probabilities
heatmap(project_2d(triplets, periodic=true),
        ticks=:auto, colorbartitle="Triplet probability",
        color=cgrad(:gist_heat, rev=true), clims=(0,0.0035), scale=:log10)
png(out_dir*"measured_probabilities_2D_average")


for (index,prediction) in enumerate(["ideal","loop_extr","pairwise", "quad_hamiltonian"])
        file_prefix=out_dir*prediction
        if prediction== "ideal"
                pred_triplets=ideal(P)
        elseif prediction=="quad_hamiltonian"
                pred_triplets=ideal(P)
                pred_triplets.*=f_factors
                pred_triplets[pred_triplets.>1].=1
        elseif prediction=="loop_extr"
                pred_triplets=loop_extr(P)
        elseif prediction=="pairwise"
                pred_triplets=pairwise_int(P,P0)
                #save the effective energy map
                heatmap(-log.(P./P0), clims=(-5,5),
                        color=cgrad(:coolwarm),
                        colorbartitle="Pairwise energy estimate")
                png(file_prefix*"_energy_estimate")
                display(current())
        else
                error("Unknown prediction label $prediction")
        end

        #Save the P_3(s) curve for later
        P_3_s[:,index+1]=triplets_1d(pred_triplets)

        #Heatmaps below. Color set for each plot type
        heatmap(pred_triplets[:,:,bait_point], clims=(0,0.0035),
                color=cgrad(:gist_heat, rev=true),
                colorbartitle="Triplet probability")
        png(file_prefix*"_triplets_bait_$(bait_point)")

        p_values=p_vals(triplets[:,:,bait_point],pred_triplets[:,:,bait_point],num_samp)

        heatmap(p_values, clims=(0,0.1), color=cgrad(:hot),
                colorbartitle="P-value")
        png(file_prefix*"_p_vals_bait_$(bait_point)")

        heatmap(adjust_array(p_values,BenjaminiHochberg()),
                colorbartitle="Benjamini-Hochberg adj. p-value",
                color=cgrad(:tempo,rev=true), clims=(0,0.3))
        png(file_prefix*"_BH_adj_p_bait_$(bait_point)")

        zs=z_scores(triplets,pred_triplets,num_samp)

        heatmap(zs[:,:,bait_point], colorbartitle="Z-score",
                color=cgrad(:coolwarm), clims=(-4,4))
        png(file_prefix*"_z_scores_bait_$(bait_point)")

        #2D average z-scores
        heatmap(project_2d(zs), scale=:log10,
                ticks=:auto, colorbartitle="Mean z-score",
                color=cgrad(:coolwarm), clims=(-4,4))
                png(file_prefix*"_z_scores_2D_average")

        heatmap(project_2d(pred_triplets, periodic=true),
                ticks=:auto, colorbartitle="Triplet probability",
                color=cgrad(:gist_heat, rev=true), clims=(0,0.0035), scale=:log10)
        png(file_prefix*"_probabilities_2D_average")
end

#Save all the P_3(s) curves
xs=(4*4:4:4*(size(P_3_s)[1]+3))*10 #axis in kb
plot(xs,P_3_s[:,1:4], label=["Condensin data" "Ideal polymer" "Loop-extruder" "Pairwise interaction"],
        xlabel="Genomic length largest loop (kb)", scale=:log10, ticks=:auto,
        ylabel="Mean triplet probability", size=[640,500])
png(out_dir*"P_3_s_curves")

plot(xs,P_3_s[:,[1,5]], label=["Condensin data" "Quadratic Hamiltonian" "Loop-extruder"],
        xlabel="Genomic length largest loop (kb)", scale=:log10, ticks=:auto,
        ylabel="Mean triplet probability", size=[640,500])
png(out_dir*"quad_hamiltonian_P_3_s_curves")