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

using HDF5, DelimitedFiles

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
close(fid)

writedlm(out_contacts,contacts)
