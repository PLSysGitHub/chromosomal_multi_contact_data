"""

This script generates P_0(i,j) maps required for the use of
the pariwise interaction formula, by simulating a random walk in confinement.

The output is stored in the folder "./Contact_files/",
where it can be accessed by the example script for generating
triplet predictions.

Takes maybe 10 minutes to run with high number of samples. Ready
made file is included.

"""

using DelimitedFiles
using ProgressMeter
using LinearAlgebra

#Helper function returns uniformly distr. point in ball of radius R
function random_start_site(R)
    in_bounds=false
    x=zeros(3)
    while !in_bounds
        x=[rand() for i in 1:3].*R #uniform in box
        in_bounds=norm(x)<R
    end
    return x
end

#Helper function returns next monomer site from x, in ball of radius R
function next_site_FJC(x,R)
    in_bounds=false
    δ=[0,0,0]
    while !in_bounds
        δ=randn(3)
        δ./=norm(δ)
        in_bounds=norm(x.+δ)<R
    end
    return x.+δ
end

"""
Simulate random walk for N_sims steps, return
the contact matrix, where contacts are defined for
every nth monomer being within rcl of each other
"""
function contacts_FJC(N::Int64,R::Number,N_sims::Int64, n::Int64, contact_radii::Vector{Float64})
    contacts=Dict(Pair.(contact_radii, [zeros(Int(N/n),Int(N/n)) for i in 1:length(contact_radii)]))
    pos=zeros(3,N-n) #-n because last n-1 monomers not physical

    pos[:,1]=random_start_site(R)
    @showprogress for i in 1:N_sims
        for j in 2:N-n
            pos[:,j].=next_site_FJC(pos[:,j-1],R)
            if j%n==1 #corresponds to bin
                for k in 1:n:j-n
                    for rcl in contact_radii
                        if norm(pos[:,k].-pos[:,j])<rcl
                            contacts[rcl][Int((k-1)/n+1),Int((j-1)/n+1)]+=1
                        end
                    end
                end
            end
        end
        #reset the simulation
        pos[:,1].=random_start_site(R)
    end
    for rcl in contact_radii
        contacts[rcl]./=N_sims
        contacts[rcl].+=transpose(contacts[rcl])
    end

    return contacts
end

"""
Turn contact map to P(s) curve, return contact map
where entries taken straight from P(s)
"""
function average_contact_map(M)
    N=size(M)[1]
    ps=zeros(N)
    counts=zeros(N)
    for i in 1:N
        for j in i+1:N
            ps[j-i]+=M[j,i]
            counts[j-i]+=1
        end
    end
    ps./=counts
    av_M=zeros(size(M))
    for i in 1:N
        for j in 1:N
            if j!=i
                av_M[j,i]=ps[abs(j-i)]
            end
        end
    end
    return av_M
end


contact_folder="./Contact_files/"
n=10                 #number of segments to map one data bin to
b=377.1669/sqrt(n)   #monomer size
N=65                 #total number of bins
R=1674.905 / 2. / b  #radius of confinement, simulation units
num_sims=1000000     #total number sampled
contact_radii=[100, 125, 150, 175, 200]

#Simulations with ProgressMeter
contacts=contacts_FJC(N*n, R, num_sims, n, contact_radii./b) #returns dict with contact maps for different radii

for contact_r in  contact_radii#nm
    #Generate the name for the output file
    out_file_name= contact_folder*"P0_linear_n$(n)_b$(b)_N$(N)_R$(R)_contact_r_$contact_r.txt"
    #Save
    av_contacts=average_contact_map(contacts[contact_r/b])
    writedlm(out_file_name,av_contacts)
end
