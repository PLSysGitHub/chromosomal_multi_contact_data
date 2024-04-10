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

out_contacts="Contact_files/bintu_IMR90.txt"
out_triplets="Triplet_files/bintu_IMR90.h5"
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
r=zeros(N,N,N) #correlation of distances, according to Polovnikov et.al. 2019
d_squared=zeros(N,N)
num_samples_3=zeros(N,N,N)
num_samples_2=zeros(N,N)
@showprogress 1 "Analysing samples..." for cell in 0:num_cells-1
    positions=data[cell*N+1:(cell+1)*N,:]
    ds=pairwise(Euclidean(),positions,dims=1)
    for i in 1:N
        for j in i+1:N
            #Are i,j in contact?
            if !isnan(ds[j,i])
                d_squared[j,i]+=ds[j,i]^2
                num_samples_2[j,i]+=1

                if ds[j,i]<contact_R
                    contacts[j,i]+=1
                end
    
                for k in j+1:N
                    if !any(isnan, [ds[j,k],ds[k,i]]) #if any of the distances are NaN, skip
                        num_samples_3[k,j,i]+=1
                        r[k,j,i]+=ds[j,i]*ds[k,j]
    
                        ds_triple=sort([ds[k,i],ds[j,i],ds[j,k]])
    
                        if ds_triple[1]<contact_R && ds_triple[2]<contact_R #if two are in contact, call it a triplet
                            triplets[k,j,i]+=1
                        end
                    end
                end
            end
        end
    end
end

triplets./=num_samples_3
contacts./=num_samples_2
r./=num_samples_3
d_squared./=num_samples_2

for i in 1:N
    for j in i+1:N
        for k in j+1:N
            r[k,j,i]/=sqrt(d_squared[j,i]*d_squared[k,j])
        end
    end
end
f=(1 .- r.^2).^(-3/2) #N,N,N array of factors f, see Polovnikov et.al. 2019

#Symmetrise arrays
for i in 1:N
    for j in i+1:N
        contacts[i,j]=contacts[j,i]
        num_samples_2[i,j]=num_samples_2[j,i]
        for k in j+1:N
            f[i,j,k]=f[k,j,i]
            f[i,k,j]=f[k,j,i]
            f[j,i,k]=f[k,j,i]
            f[j,k,i]=f[k,j,i]
            f[k,i,j]=f[k,j,i]
            triplets[i,j,k]=triplets[k,j,i]
            triplets[i,k,j]=triplets[k,j,i]
            triplets[j,i,k]=triplets[k,j,i]
            triplets[j,k,i]=triplets[k,j,i]
            triplets[k,i,j]=triplets[k,j,i]
            num_samples_3[i,j,k]=num_samples_3[k,j,i]
            num_samples_3[i,k,j]=num_samples_3[k,j,i]
            num_samples_3[j,i,k]=num_samples_3[k,j,i]
            num_samples_3[j,k,i]=num_samples_3[k,j,i]
            num_samples_3[k,i,j]=num_samples_3[k,j,i]
        end
    end
end

#save the data
fid=h5open(out_triplets,"w")
fid["triplets"]=triplets
fid["num_samples"]=num_cells
fid["num_samples_2"]=num_samples_2 #save separately, to account for NaN values
fid["num_samples_3"]=num_samples_3 #save separately, to account for NaN values
fid["f_factors"]=f
close(fid)

writedlm(out_contacts,contacts)