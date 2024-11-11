"""
This file contains functions for calculating contact triplet predictions as
described in the preprint
"Predicting chromosomal multi-contact data using pairwise contact frequencies"

Many functions have a variable periodic=true, which can should be set to false
when not considering bacterial chromosomes.
"""


using HypothesisTests, MultipleTesting, StatsBase, DelimitedFiles, HDF5

"""
Helper function gives genomic length of loop (a,b) on
a circular chromosome of length N
"""
d_periodic(a,b,N)=min(abs(a-b),N-abs(a-b))

"""
Given three points, return sorted triplet (i,j,k) so that
(i,k) corresponds to the largest loop
"""
function sort_periodic(a,b,c, N)
    ab=d_periodic(a, b, N)
    bc=d_periodic(b, c, N)
    ac=d_periodic(a, c, N)

    maxd=maximum([ab,bc,ac])

    if ab==maxd
        return [a,c,b]
    elseif bc==maxd
        return [b,a,c]
    else
        return [a,b,c]
    end
end

"""
Given three points, return the lengths of the two smaller loops
x,y and s=x+y
"""
function loop_lengths_periodic(a,b,c,N)
    ab=d_periodic(a, b, N)
    bc=d_periodic(b, c, N)
    ac=d_periodic(a, c, N)

    dists=sort([ab,bc,ac])
    dists[3]=dists[1]+dists[2]
    return dists
end

"""
Given pair-wise contact frequencies P(i,j)
Return 3D prediction for contact triplets

    P(i,j,k) ≈ P(i,j)P(j,k)+P(j,k)P(k,i)+P(k,i)P(i,j)
"""
function loop_extr_first_order(contacts;periodic=true)
    N=size(contacts)[1]
    triplets=zeros(N,N,N)
    for i in 1:N
        for j in 1:N
            for k in 1:N
                if periodic
                    if d_periodic(i,j,N)>1 && d_periodic(j,k,N)>1 && d_periodic(i,k,N)>1
                        p=min(1,contacts[i,j]*contacts[j,k]+contacts[j,k]*contacts[k,i]+contacts[k,i]*contacts[i,j])
                        triplets[k,j,i]=p
                    end
                else
                    if i!=j && k!=i && j!=k
                        p=min(1,contacts[i,j]*contacts[j,k]+contacts[j,k]*contacts[k,i]+contacts[k,i]*contacts[i,j])
                        triplets[k,j,i]=p
                    end
                end
            end
        end
    end
    return triplets
end

"""
Given pair-wise contact frequencies P(i,j)
and the non-interacting background P_0(i,j)
Return 3D prediction for contact triplets

    P(i,j,k) ≈ P(i,j)P(j,k)+
    P(j,k)P(k,i)+
    P(k,i)P(i,j)-
    2P(i,j)P(j,k)P(k,i)

"""
function loop_extr(contacts;periodic=true)
    N=size(contacts)[1]
    triplets=zeros(N,N,N)
    for i in 1:N
        for j in 1:N
            for k in 1:N
                if periodic
                    if d_periodic(i,j,N)>1 && d_periodic(j,k,N)>1 && d_periodic(i,k,N)>1
                        p=min(1,contacts[i,j]*contacts[j,k]+contacts[i,k]*contacts[k,j]+contacts[j,i]*contacts[i,k]-2*contacts[i,j]*contacts[j,k]*contacts[k,i])
                        triplets[k,j,i]=p
                    end
                else
                    if i!=j && k!=i && j!=k
                        p=min(1,contacts[i,j]*contacts[j,k]+contacts[i,k]*contacts[k,j]+contacts[j,i]*contacts[i,k]-2*contacts[i,j]*contacts[j,k]*contacts[k,i])
                        triplets[k,j,i]=p
                    end
                end
            end
        end
    end
    return triplets
end

"""
Given pair-wise contact frequencies P(i,j)
Return prediction for contact triplets.
For i<j<k, the ideal polymer prediction is:

    P(i,j,k) ≈ P(i,j)P(j,k)
"""
function ideal(contacts;periodic=true)
    N=size(contacts)[1]
    triplets=zeros(N,N,N)
    for i in 1:N
        for j in 1:N
            for k in 1:N
                if periodic
                    if d_periodic(i,j,N)>1 && d_periodic(j,k,N)>1 && d_periodic(k,i,N)>1
                        a,b,c=sort_periodic(i,j,k,N)
                        p=contacts[a,b]*contacts[b,c]
                        triplets[k,j,i]=p
                    end
                else
                    if i!=j && k!=i && j!=k
                        a,b,c=sort([i,j,k])
                        p=contacts[a,b]*contacts[b,c]
                        triplets[k,j,i]=p
                    end
                end
            end
        end
    end
    return triplets
end

"""
Given pair-wise contact frequencies P(i,j)
and corresponding probabilities P_0(i,j)
Return 3D prediction for contact triplets

    P(i,j,k) ≈ P(i,j)P(j,k)P(k,i)/P_0(k,i)

"""
function pairwise_int(contacts, non_interacting;periodic=true)
    @assert size(contacts)==size(non_interacting) "Different sized arrays!"

    N=size(contacts)[1]
    triplets=zeros(N,N,N)
    for i in 1:N
        for j in 1:N
            for k in 1:N
                if periodic
                    if d_periodic(i,j,N)>1 && d_periodic(j,k,N)>1 && d_periodic(k,i,N)>1
                        a,b,c=sort_periodic(i,j,k,N)
                        p=min(1,contacts[a,b]*contacts[b,c]*contacts[c,a]/non_interacting[c,a])
                        triplets[k,j,i]=p
                    end
                else
                    if i!=j && k!=i && j!=k
                        a,b,c=sort([i,j,k])
                        p=min(1,contacts[a,b]*contacts[b,c]*contacts[a,c]/non_interacting[a,c])
                        triplets[k,j,i]=p
                    end
                end
            end
        end
    end
    triplets[isinf.(triplets)].=NaN
    return triplets
end


"""
Given 3D array of triplet frequencies P(i,j,k)
Return 1D projection P(s), where s=x+y, and x
and y are the lengths of the two smallest loops

Bins with zero counts are output as NaN
"""
function triplets_1d(triplets::Array{Float64,3};periodic=true)
    N=size(triplets)[1]
    if periodic
        result=zeros(floor(Int64,2*N/3)) #largest triplet when all points evenly spaced
        counts=zeros(floor(Int64,2*N/3))
        for i in 1:N
            for j in i+2:N
                for k in j+2:N
                    x,y,s=loop_lengths_periodic(i,j,k,N) #s=x+y
                    result[s]+=triplets[k,j,i]
                    counts[s]+=1
                end
            end
        end
        result./=counts
        result[result.==Inf].=NaN
        result[iszero.(result)].=NaN #for log plotting

        #shortest triplet is size 4.
        #Upper bound because for larger s doesn't correspond to genomic length
        return result[4:floor(Int64,N/2)]
    else
        result=zeros(N-1)
        counts=zeros(N-1)
        for i in 1:N
            for j in i+1:N
                for k in j+1:N
                    s=k-i
                    if !isnan(triplets[k,j,i])
                        result[s]+=triplets[k,j,i]
                        counts[s]+=1
                    end
                end
            end
        end
        result./=counts
        result[result.==Inf].=NaN
        result[iszero.(result)].=NaN #for log plots

        return result[2:end]
    end
end

function P_s(contacts;periodic=true)
    N=size(contacts)[1]
    P_s=zeros(N-1)
    counts=zeros(N-1)
    if !periodic
        for i in 1:N
            for j in i+1:N
                s=j-i
                P_s[s]+=contacts[i,j]
                counts[s]+=1
            end
        end
    else
        for i in 1:N
            for j in i+1:N
                s=d_periodic(i,j,N)
                P_s[s]+=contacts[i,j]
                counts[s]+=1
            end
        end
    end
    return P_s./counts
end

"""
Given 3D array M(i,j,k)
Return 2D projection A(i,j), where
x and y are the genomic lengths of the
two smallest loops in (i,j,k)

Bins with zero counts are output as NaN
"""
function project_2d(M::Array{Float64,3};periodic=true)
    N=size(M)[1]
    if periodic
        #max size of smaller loop is for (1,2,N/2)
        A=zeros(floor(Int64,N/3),floor(Int64,N/2))
        counts=zeros(floor(Int64,N/3),floor(Int64,N/2))
        for i in 1:N
            for j in i+2:N
                for k in j+2:N
                    x,y,s=loop_lengths_periodic(i,j,k,N)
                    if x>1 && y>1
                        A[x,y]+=M[k,j,i]
                        counts[x,y]+=1
                    end
                end
            end
        end
        A[iszero.(A)].=NaN
        A./=counts
        return A[3:end,3:end]
    else
        A=zeros(floor(Int64,N/2),N-1)
        counts=zeros(floor(Int64,N/2),N-1)
        for i in 1:N
            for j in i+1:N
                for k in j+1:N
                    y,x=sort([j-i,k-j])
                    if !isnan(M[k,j,i])
                        A[y,x]+=M[k,j,i]
                        counts[y,x]+=1
                    end
                end
            end
        end

        A[iszero.(counts)].=NaN
        A[iszero.(A)].=NaN
        return A./counts
    end
end

"""
Calculate elementwise z-scores presuming
a binomial distribution with probabilities
given in pred.

Data should be frequency data.
"""
function z_scores(data,pred,N_samples)
    zs=sqrt(N_samples)*(data.-pred)./sqrt.(pred.*(1 .- pred)) #difference / std
    zs[isinf.(zs)].=NaN # if prediction was 0
    return zs
end

"""
Calculate elementwise two-tailed p-values presuming
a binomial distribution with probabilities
given in pred.

Data should be frequency data.
"""
function p_vals(data,pred,N_samples::Int64)
    ps=zeros(size(data))
    non_nan=.!isnan.(pred)
    ps[non_nan].=pvalue.(BinomialTest.(round.(Int64,N_samples*data[non_nan]),N_samples,pred[non_nan]))
    ps[isnan.(data)].=NaN
    return ps
end

"""
MultipleTesting adjustment of an array,
rather than vector, of p-values
"""
function adjust_array(pValues,method)
    original_size=size(pValues)
    return reshape(adjust(vec(pValues),method),original_size)
end

"""
Function that returns matrix with A in lower half
and B in upper half. For plotting two matrices at once.
"""
function half_half(A,B)
    N=size(A)[1]
    halves=zeros(N,N)
    for i in 1:N
        for j in i+1:N
            halves[i,j]=A[i,j]
            halves[j,i]=B[i,j]
        end
    end
    return halves
end
