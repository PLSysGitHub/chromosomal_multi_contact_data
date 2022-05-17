"""
Simple script to extract triplet frequencies from
position data. Can be used to extract triplet frequencies
from Bintu et.al. files.
"""

using DelimitedFiles #Bintu et.al. .txt files
using HDF5 #For saving triplet array
using Distances #Effective distance matrix
using LinearAlgebra #norm
#=
File format:
cell, locus, x, y, z, (cell cycle phase)
=#

out_contacts="Contact_files/bintu_IMR90_test.txt"
out_triplets="Triplet_files/bintu_IMR90_test.h5"
in_file_name="Extract_triplets/IMR90_chr21-28-30Mb_cell cycle.csv"

if !isfile(in_file_name)
    println("Downloading data file...")
    using Downloads
    Downloads.download("https://github.com/BogdanBintu/ChromatinImaging/raw/master/Data/IMR90_chr21-28-30Mb_cell%20cycle.csv",in_file_name)
end

contact_R=150 #nm
data=readdlm(in_file_name,',', skipstart=2)
N=Int(maximum(data[end,2]))
num_cells=Int(length(data[:,1])/N)
data=convert.(Float64, data[:,3:5]) #positions for monomers in each cell

triplets=zeros(N,N,N)
contacts=zeros(N,N)

@time for cell in 0:num_cells-1
    positions=data[cell*N+1:(cell+1)*N,:]
    ds=pairwise(Euclidean(),positions,dims=1)
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
                    ds_triple=sort([ds[k,i],ds[j,i],ds[j,k]])

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
end

triplets./=num_cells
contacts./=num_cells

#save the data
fid=h5open(out_triplets,"w")
fid["triplets"]=triplets
close(fid)

writedlm(out_contacts,contacts)
