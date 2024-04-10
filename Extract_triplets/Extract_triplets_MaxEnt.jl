"""
The MaxEnt simulations are for a lattice polymer.
The simulations include a built-in XXX array, which lists
all the monomers [i,j,k,l,...] that are in contact at any
time.

The easiest way to save contact information is to simply output
XXX to a file the required amount of times. This allows for saving
far more samples than if saving full 3D configurations.

This directory includes such a file, "MaxEnt_contact_data.txt"

The following script uses the file to extract both contact and contact
triplet frequencies. Higher order contacts and other contact structure
frequencies can also be found efficiently.

"""


using HDF5, DelimitedFiles, ProgressMeter, Distances, LinearAlgebra

"""
Return the non-normalised contact map calculated from
contact list file

File should have rows
i j k n
m l

etc, where each row corresponds to contacts between all
monomer indices on it. Empty row represents end of one
configuration.

n is the number to divide indices by (C++ MaxEnt simulations
use n=4)
"""
function sample_contacts(file_name::String,pol_length, n)
    sampled_hi_c=zeros(pol_length,pol_length)
    num_samples=0
    open(file_name) do file
        for ln in eachline(file)
            cluster=parse.(Int,split(ln))
            cluster.=((cluster)./n .+ 1) #C++ indexing starts from 0, unlike Julia
            if length(cluster)==0
                num_samples+=1
            else
                for i in 1:length(cluster)
                    for j in i+1:length(cluster)
                        sampled_hi_c[cluster[i],cluster[j]]+=1
                    end
                end
            end
        end
    end

    # Set self- and neighbour-contacts to zero
    for i in 1:pol_length
        sampled_hi_c[i%pol_length+1,i]=0
        sampled_hi_c[i,i]=0
        sampled_hi_c[i,i%pol_length+1]=0
    end

    sampled_hi_c.+=transpose(sampled_hi_c)
    sampled_hi_c./=num_samples

    return sampled_hi_c, num_samples
end

"""
Returns a 3D array of contact triplet frequencies P(i,j,k),
using same contact list file as sample_contacts.

n is the number to divide indices by (C++ MaxEnt simulations
use n=4)
"""
function sample_triplets(file_name::String, pol_length, n)
    density_3d=zeros(pol_length,pol_length, pol_length)
    num_samples=0
    open(file_name) do file
        for ln in eachline(file)
            cluster=parse.(Int,split(ln))
            cluster.=((cluster)./n .+ 1)

            if length(cluster)==0
                num_samples+=1
            elseif length(cluster)>2
                for i in cluster
                    for j in cluster
                        if i!=j && i!=j%pol_length+1 && j!=i%pol_length+1
                            for k in cluster
                                if j!=k && k!=i && k!=j%pol_length+1 && j!=k%pol_length+1 && i!=k%pol_length+1 && k!=i%pol_length+1
                                    density_3d[i,j,k]+=1
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return density_3d./num_samples
end

#adapt above for quadruplets etc


N=405 #Number of bins
n=4 #number of monomers a bin is mapped to. 4 by default in MaxEnt program

out_triplets="Triplet_files/MaxEnt_triplets.h5"
out_contacts="Contact_files/c_crescentus_hic.txt"
contact_file="Extract_triplets/MaxEnt_contact_data.txt"

#=
Note; C++ simulations can be used to give contact map directly
This is what was done for the non-interacting polymer
=#

contacts, num_samples=sample_contacts(contact_file, N,n)
triplets=sample_triplets(contact_file, N,n)

#save the data
fid=h5open(out_triplets,"w")
fid["triplets"]=triplets
fid["num_samples"]=num_samples
close(fid)

writedlm(out_contacts,contacts)

#For Polovnikov et.al. 2019, we need distance correlations
r=zeros(N,N,N) #correlation of distances, according to Polovnikov et.al. 2019
d_squared=zeros(N,N)
f_name="Triplet_files/MaxEnt_configs_3D.h5"
num_3d_samples=0
h5open(f_name,"r") do fid
    ks=filter(x->isnumeric(x[1]), keys(fid))
    global num_3d_samples=length(ks)
    @showprogress 1 "Looping through positions to calculate f-factors..." for i in ks
        data=read(fid,"$i") 
        ds=pairwise(Euclidean(),data,dims=2)
        d_squared.+=ds.^2
        r .+= reshape(ds, N, N, 1) .* reshape(ds, 1, N, N) #dij*djk
    end
end

r./=num_3d_samples #now <r_ij r_jk> for all i,j,k
d_squared./=num_3d_samples #now <r_ij^2> for all i,j
#Divide to get r = <r_ij r_jk> / sqrt(<r_ij^2> <r_jk^2>)
r./= sqrt.(reshape(d_squared, N, N, 1) .* reshape(d_squared, 1, N, N))
f=(1 .- r.^2).^(-3/2) #N,N,N array of factors f, see Polovnikov et.al. 2019

#symmetrise for i<j<k
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
fid=h5open(out_triplets,"r+")
fid["f_factors"]=f
close(fid)