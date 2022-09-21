"""
Simple script to extract triplet frequencies from
position data. Can be used to extract triplet frequencies
from loop extruder simulated files from script of Brandão et.al.

Takes a while to run.
"""

using HDF5 #Brandão et.al. .h5 files
using DelimitedFiles #For saving contact map
using LinearAlgebra #Norm
using Distances

#Folder with the .h5 files from the loop extruder simulations
folder="Extract_triplets/Example_Brandao/"
out_triplets="Triplet_files/condensin_triplets.h5"
out_contacts="Contact_files/condensin_hic.txt"

## Similar for non-interacting. Note no need to save triplets!
# folder="Extract_triplets/Example_Brandao_non_int/"
# out_contacts="Contact_files/P0_loop_extruder.txt"


file_names=readdir(folder, join=true)
file_names=file_names[contains.(file_names,"blocks_")] #don't want any other files

N=404
contacts=zeros(N,N)
triplets=zeros(N,N,N)
num_samples=0
contact_R=5

for f_name in file_names
    fid=h5open(f_name, "r")

    for cell in keys(fid)
        data=read(fid,"$cell/pos")[:,1:10:end] #for 10 kb bins use every 10th
        ds=pairwise(Euclidean(),data,dims=2)
        for i in 1:N
            for j in i+1:N
                #Are i,j in contact?
                if ds[j,i]<contact_R
                    contacts[j,i]+=1
                    contacts[i,j]+=1
                end

                #Only look for triplets if can occur
                if ds[j,i]<2*contact_R
                    for k in j+1:N
                        ds_triple=sort([ds[k,j],ds[j,i],ds[k,i]])

                        if ds_triple[1]<contact_R && ds_triple[2]<contact_R #if two are in contact, call it a triplet
                            triplets[j,k,i]+=1
                            triplets[k,j,i]+=1
                            triplets[k,i,j]+=1
                            triplets[i,k,j]+=1
                            triplets[j,i,k]+=1
                            triplets[i,j,k]+=1
                        end

                    end
                end
            end
        end
        num_samples+=1
    end
end
contacts./=num_samples
triplets./=num_samples

#save the data
fid=h5open(out_triplets,"w")
fid["triplets"]=triplets
close(fid)

writedlm(out_contacts,contacts)
