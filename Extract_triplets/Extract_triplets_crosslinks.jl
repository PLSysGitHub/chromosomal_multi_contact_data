using JLD2, DelimitedFiles, HDF5


"""
The crosslinker simulations are for a lattice polymer.
The simulations include a built-in "contacts" array, which lists
all the monomers [i,j,k,l,...] that are in contact at any
time.

The easiest way to save contact information is to simply output
"contacts" to a file the required amount of times. This allows for saving
far more samples than if saving full 3D configurations.

The directory "Crosslinks" includes such txt files for different crosslinker
potentials. The files have the name "mu_[mu]_[random/cos]_[A]", where mu is the
chemical potential for crosslinkers, random/cos describes whether the binding energies
are sampled uniformly from [-A,A] or given by A*cos((i))

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

n is the number to divide indices by
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

n is the number to divide indices by
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
n=1 #number of monomers a bin is mapped to.

#=
Note; C++ simulations can be used to give contact map directly
This is what was done for the non-interacting polymer
=#

num_samples=Dict{String, Int}()

files=readdir("Extract_triplets/Crosslinks")
files=files[contains.(files,".txt")]

for f in files
    contact_file="Extract_triplets/Crosslinks/$(f)"

    contacts, num_s=sample_contacts(contact_file, N,n)

    num_samples[f[1:end-4]]=num_s
    writedlm("Contact_files/Crosslinks/$f",contacts)

    triplets=sample_triplets(contact_file, N,n)
    fid=h5open("Triplet_files/Crosslinks/$(f[1:end-4]).h5","w")
    fid["triplets"]=triplets
    close(fid)
end

jldsave("./Extract_triplets/Crosslinks/sample_counts.jld2"; num_samples)
