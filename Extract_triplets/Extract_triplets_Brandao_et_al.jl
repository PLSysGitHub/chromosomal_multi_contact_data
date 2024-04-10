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
# out_contacts="Contact_files/P0_condensin.txt"


file_names=readdir(folder, join=true)
file_names=file_names[contains.(file_names,"blocks_")] #don't want any other files

N=404
contacts=zeros(N,N)
triplets=zeros(N,N,N)
r=zeros(N,N,N) #correlation of distances, according to Polovnikov et.al. 2019
d_squared=zeros(N,N)
num_samples=0
contact_R=5

@showprogress 1 "Looping through files" for file_name in file_names
    fid=h5open(file_name, "r")
    for cell in keys(fid)
        data=read(fid,"$cell/pos")[:,1:10:end] #for 10 kb bins use every 10th
        ds=pairwise(Euclidean(),data,dims=2)
        d_squared.+=ds.^2
        r .+= reshape(ds, size(ds, 1), size(ds, 2), 1) .* reshape(ds, 1, size(ds, 2), size(ds, 1))
        #Check for contacts and triplets
        for i in 1:N
            for j in i+1:N
                #Are i,j close enough for a triplet?
                if ds[j,i]<=2*contact_R
                    #Are i,j in contact?
                    if ds[j,i]<contact_R
                        contacts[j,i]+=1
                        contacts[i,j]+=1
                    end

                    for k in j+1:N
                        ds_triple=sort([ds[k,j],ds[j,i],ds[k,i]])
                        if ds_triple[1]<contact_R && ds_triple[2]<contact_R #if two are in contact, call it a triplet
                            triplets[k,j,i]+=1
                        end
                    end
                end
            end
        end
        global num_samples+=1
    end
    close(fid)
end
contacts./=num_samples
triplets./=num_samples
for i in 1:N
    for j in i+1:N
        for k in j+1:N
            triplets[j,k,i]=triplets[k,j,i]
            triplets[k,i,j]=triplets[k,j,i]
            triplets[i,j,k]=triplets[k,j,i]
            triplets[i,k,j]=triplets[k,j,i]
            triplets[j,i,k]=triplets[k,j,i]
        end
    end
end

r./=num_samples #now <r_ij r_jk> for all i,j,k
d_squared./=num_samples #now <r_ij^2> for all i,j
#Divide to get r = <r_ij r_jk> / sqrt(<r_ij^2> <r_jk^2>)
r./= sqrt.(reshape(d_squared, size(d_squared, 1), size(d_squared, 2), 1) .* reshape(d_squared, 1, size(d_squared, 2), size(d_squared, 1)))
f=(1 .- r.^2).^(-3/2) #N,N,N array of factors f, see Polovnikov et.al. 2019

for i in 1:N
    for j in i+1:N
        for k in j+1:N
            f[k,j,i]=f[i,j,k]
            f[i,k,j]=f[i,j,k]
            f[j,i,k]=f[i,j,k]
            f[j,k,i]=f[i,j,k]
            f[k,i,j]=f[i,j,k]
        end
    end
end

#save the data
out_fid=h5open(out_triplets,"w")
out_fid["triplets"]=triplets
out_fid["f_factors"]=f
out_fid["num_samples"]=num_samples
close(out_fid)

writedlm(out_contacts,contacts)
