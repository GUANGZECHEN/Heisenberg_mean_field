using LinearAlgebra
using DelimitedFiles
using Statistics
#using Plots
using DelimitedFiles
include("Geometry.jl")
include("pi-flux.jl")
include("rearrange_results.jl")
include("matrix_initialization.jl")

function mag_H(J1,J2,J3,M,N,R,inter_vector) # M:magnetization matrix, N:total number of sites
    H=zeros(Complex,2*N,2*N)
    n_inter=size(inter_vector,1)
    for nn=1:n_inter
        r=[]
        for ii=1:N
            push!(r,inter_vector[nn])
        end
        R2=R+r
    
        for i=1:N
            for j=1:N
                if nearest(R[i],R2[j])   # magnetization from nearest J1 interaction
                    H[2*i-1,2*i-1]=H[2*i-1,2*i-1]+1/4*J1*M[2*j-1,2*j-1]
                    H[2*i,2*i]=H[2*i,2*i]+1/4*J1*M[2*j,2*j]
                    H[2*i,2*i-1]=H[2*i,2*i-1]+1/4*J1*M[2*j-1,2*j]
                    H[2*i-1,2*i]=H[2*i-1,2*i]+1/4*J1*conj(M[2*j-1,2*j])

                    H[2*j-1,2*j-1]=H[2*j-1,2*j-1]+1/4*J1*M[2*i-1,2*i-1]
                    H[2*j,2*j]=H[2*j,2*j]+1/4*J1*M[2*i,2*i]
                    H[2*j,2*j-1]=H[2*j,2*j-1]+1/4*J1*M[2*i-1,2*i]
                    H[2*j-1,2*j]=H[2*j-1,2*j]+1/4*J1*conj(M[2*i-1,2*i])
                elseif next_nearest(R[i],R2[j])   # magnetization from next_nearest J2 interaction
                    H[2*i-1,2*i-1]=H[2*i-1,2*i-1]+1/4*J2*M[2*j-1,2*j-1]
                    H[2*i,2*i]=H[2*i,2*i]+1/4*J2*M[2*j,2*j]
                    H[2*i,2*i-1]=H[2*i,2*i-1]+1/4*J2*M[2*j-1,2*j]
                    H[2*i-1,2*i]=H[2*i-1,2*i]+1/4*J2*conj(M[2*j-1,2*j])

                    H[2*j-1,2*j-1]=H[2*j-1,2*j-1]+1/4*J2*M[2*i-1,2*i-1]
                    H[2*j,2*j]=H[2*j,2*j]+1/4*J2*M[2*i,2*i]
                    H[2*j,2*j-1]=H[2*j,2*j-1]+1/4*J2*M[2*i-1,2*i]
                    H[2*j-1,2*j]=H[2*j-1,2*j]+1/4*J2*conj(M[2*i-1,2*i])
                elseif next_next_nearest(R[i],R2[j])   # magnetization from next_next_nearest J3 interaction
                    H[2*i-1,2*i-1]=H[2*i-1,2*i-1]+1/4*J3*M[2*j-1,2*j-1]
                    H[2*i,2*i]=H[2*i,2*i]+1/4*J3*M[2*j,2*j]
                    H[2*i,2*i-1]=H[2*i,2*i-1]+1/4*J3*M[2*j-1,2*j]
                    H[2*i-1,2*i]=H[2*i-1,2*i]+1/4*J3*conj(M[2*j-1,2*j])

                    H[2*j-1,2*j-1]=H[2*j-1,2*j-1]+1/4*J3*M[2*i-1,2*i-1]
                    H[2*j,2*j]=H[2*j,2*j]+1/4*J3*M[2*i,2*i]
                    H[2*j,2*j-1]=H[2*j,2*j-1]+1/4*J3*M[2*i-1,2*i]
                    H[2*j-1,2*j]=H[2*j-1,2*j]+1/4*J3*conj(M[2*i-1,2*i])
                end
            end
        end
    end
    
    if norm(H-adjoint(H))<0.0001
        return (H+adjoint(H))/2
    else
        println("mag_H is not Hermitian")
    end
end

function hop_H(J1,J2,J3,Chi,N,R,inter_vector) #Chi: hopping matrix
    H=zeros(Complex,2*N,2*N)
    n_inter=size(inter_vector,1)

    for nn=1:n_inter
        r=[]
        for ii=1:N
            push!(r,inter_vector[nn])
        end
        R2=R+r

        for i=1:N
            for j=1:N
                if nearest(R[i],R2[j])   # hopping from nearest J1 interaction
                    H[2*j-1,2*i-1]=H[2*j-1,2*i-1]-1/4*J1*Chi[i,j]
                    H[2*j,2*i]=H[2*j,2*i]-1/4*J1*Chi[i,j]    

                    H[2*i-1,2*j-1]=H[2*i-1,2*j-1]-1/4*J1*Chi[j,i]
                    H[2*i,2*j]=H[2*i,2*j]-1/4*J1*Chi[j,i]               
                elseif next_nearest(R[i],R2[j])   # hopping from next_nearest J2 interaction
                    H[2*j-1,2*i-1]=H[2*j-1,2*i-1]-1/4*J2*Chi[i,j]
                    H[2*j,2*i]=H[2*j,2*i]-1/4*J2*Chi[i,j]
                    
                    H[2*i-1,2*j-1]=H[2*i-1,2*j-1]-1/4*J2*Chi[j,i]
                    H[2*i,2*j]=H[2*i,2*j]-1/4*J2*Chi[j,i] 
                elseif next_next_nearest(R[i],R2[j])   # hopping from next_next_nearest J3 interaction
                    H[2*j-1,2*i-1]=H[2*j-1,2*i-1]-1/4*J3*Chi[i,j]
                    H[2*j,2*i]=H[2*j,2*i]-1/4*J3*Chi[i,j]

                    H[2*i-1,2*j-1]=H[2*i-1,2*j-1]-1/4*J3*Chi[j,i]
                    H[2*i,2*j]=H[2*i,2*j]-1/4*J3*Chi[j,i] 
                end
            end
        end
    end
    
    if norm(H-adjoint(H))<0.0001
        return (H+adjoint(H))/2
    else
        println("hop_H is not Hermitian")
    end
end

function chemical_potential(N,Mu) #Mu:chemical potential vector
    H=zeros(Complex,2*N,2*N)
    for i=1:N
        H[2*i-1,2*i-1]=Mu[i]
        H[2*i,2*i]=Mu[i]
    end
    if norm(H-adjoint(H))<0.0001
        return (H+adjoint(H))/2
    else
        println("chemical_H is not Hermitian")
    end
end

function get_H_MF(J1,J2,J3,N,M,Chi,Mu,R,inter_vector)
    H1=mag_H(J1,J2,J3,M,N,R,inter_vector)
    H2=hop_H(J1,J2,J3,Chi,N,R,inter_vector)
    #H3=chemical_potential(N,Mu)
    H=H1+H2
    return H
end

function get_Syst(J1,J2,J3,N,M,Chi,Mu,R,inter_vector,lambda,lambda_2)
    Eigvals=Float64[]
    Eigvecs=[]
    
    M0=Matrix{Complex}(I,2*N,2*N)/2

    H=lambda*get_H_MF(J1,J2,J3,N,M,Chi,Mu,R,inter_vector)+lambda_2*1/4*abs(J1)*get_H_pi_flux(R,inter_vector,N,0.4)+0*lambda_2*mag_H(J1,J2,J3,M0,N,R,inter_vector)+chemical_potential(N,Mu)
    if norm(H-adjoint(H))<0.0001
        H=(H+adjoint(H))/2
    else
        println("H_total is not Hermitian")
    end
    F=eigen(H)
           
    for l=1:2*N
        push!(Eigvals,F.values[l])
        push!(Eigvecs,F.vectors[:,l])
    end           

    return Eigvals, Eigvecs
end


