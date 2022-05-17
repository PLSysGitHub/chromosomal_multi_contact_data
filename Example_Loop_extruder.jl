"""
Example script for generating plots for predictions of contact triplets
on a bacterial chromosome.

Script calculates independent link, independent loop, and pairwise interaction predictions
Compares each to data, and saves plots in given directory

"""

include("ContactTripletPredictions.jl")
using Plots
gr(label="",size=(540,500),xticks=(0:100:400, 0:1000:4000),
        yticks=(100:100:400,1000:1000:4000), #units in kb
        foreground_color_legend = nothing, background_color_legend=nothing)

#Change names of files and directories to analyse different data
num_samp=3060
contact_file="./Contact_files/loop_extruder_hic_test.txt"
P0_file="./Contact_files/P0_loop_extruder.txt"
triplet_file="./Triplet_files/loop_extruder_triplets_test.h5" #triplet file should contain 3D triplet frequency array as "triplets"
out_dir="./Output/loop_extruder_example/"
if !isdir(out_dir) mkdir(out_dir) end

#Read in the data
P=readdlm(contact_file)
P0=readdlm(P0_file)
triplets=read(h5open(triplet_file,"r"),"triplets")

#Start pipeline
P_3_s = triplets_1d(triplets)
P_3_s = cat(P_3_s, zeros(length(P_3_s),3), dims=2) #space for predicted P_3(s) curves

#Point around which to build plots
bait_point=300

for (index,prediction) in enumerate(["ind_loop","ind_link","pairwise"])
        file_prefix=out_dir*prediction
        if prediction== "ind_loop"
                pred_triplets=ind_loop(P)
        elseif prediction=="ind_link"
                pred_triplets=ind_link(P)
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
        heatmap(project_2d(zs), scale=:log10,
                ticks=:auto, colorbartitle="Mean z-score",
                color=cgrad(:coolwarm), clims=(-4,4))
                png(file_prefix*"_z_scores_2D_average")
        display(current())
end

#Save all the P_3(s) curves
xs=(4*4:4:4*(size(P_3_s)[1]+3))*10 #axis in kb
plot(xs,P_3_s, label=["Data" "Ind. link" "Ind. loop" "Pairwise"],
        xlabel="Genomic length largest loop (kb)", scale=:log10, ticks=:auto,
        ylabel="Mean triplet probability", size=[500,400])
png(out_dir*"P_3_s_curves")
