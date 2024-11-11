"""
Example script for generating plots for predictions of contact triplets
on a bacterial chromosome with condensins.

Script calculates ideal polymer, loop-extruder, and pairwise interaction predictions
Compares each to data, and saves plots in given directory

"""

include("ContactTripletPredictions.jl")
include("Plot_functions_bacterial.jl")

#Change names of files and directories to analyse different data
num_samp=3060

for contact_r in 3:6
        contact_file="./Contact_files/condensin_hic_contact_r_$contact_r.txt"
        P0_file="./Contact_files/P0_condensin_contact_r_$contact_r.txt"
        triplet_file="./Triplet_files/condensin_triplets_contact_r_$contact_r.h5" #triplet file should contain 3D triplet frequency array as "triplets"
        out_dir="./Output/Condensin_example/Contact_radius_$contact_r/"
        mkpath(out_dir)

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
        P_3_s = cat(P_3_s, zeros(length(P_3_s),4), dims=2) #space for predicted P_3(s) curves

        #Point around which to build plots
        bait_point=200

        #heatmap(P, clims=(0,0.25), color=cgrad(:dense), colorbartitle="Contact probability")
        hic_plot(P,out_dir*"contact_map")

        #Make comparison plots for scalings
        P_S=P_s(P, periodic=true)
        plot((2:N)*10,P_S, scale=:log10, xticks=:auto, yticks=:auto, ylabel="Condensin simulation P(s)", xlabel="s")
        png(out_dir*"P_s")

        #Actual triplet probabilities
        triplet_plot(triplets, out_dir*"measured_triplet_probabilities", bait_point, label="Measured triplet probability")

        #2D average probabilities
        average_2D_plot(triplets, out_dir*"measured_triplet_probabilities_2D_average")

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
                        pred_triplets[pred_triplets.>1].=1
                        #save the effective energy map
                        heatmap(-log.(P./P0), clims=(-5,5),
                                color=cgrad(:coolwarm),
                                colorbartitle="\nPairwise energy estimate")
                        png(file_prefix*"_energy_estimate")
                else
                        error("Unknown prediction label $prediction")
                end

                #Save the P_3(s) curve for later
                P_3_s[:,index+1]=triplets_1d(pred_triplets)

                #Heatmaps below. Color set for each plot type
                triplet_plot(pred_triplets, file_prefix*"_triplets_bait_$(bait_point)", bait_point)

                p_values=p_vals(triplets[:,:,bait_point],pred_triplets[:,:,bait_point],num_samp)

                p_value_plot(p_values, file_prefix*"_p_vals_bait_$(bait_point)")

                bh_p_value_plot(p_values, file_prefix*"_BH_adj_p_bait_$(bait_point)")

                zs=z_scores(triplets,pred_triplets,num_samp)
                z_score_plot(zs, bait_point, file_prefix*"_z_scores_bait_$(bait_point)")

                #2D average z-scores
                average_2D_z_score_plot(zs, file_prefix*"_z_scores_2D_average")

                average_2D_plot(pred_triplets, file_prefix*"_probabilities_2D_average")
        end

        #Save all the P_3(s) curves
        plot_p_3_s_curves(P_3_s[:,1:4], ["Condensin sim." "Ideal" "Loop-extruder" "Pairwise"], out_dir*"P_3_s_curves")

        plot_p_3_s_curves(P_3_s[:,[1,5]], ["Condensin sim." "Quadratic Hamiltonian"], out_dir*"quad_hamiltonian_P_3_s_curves")
end