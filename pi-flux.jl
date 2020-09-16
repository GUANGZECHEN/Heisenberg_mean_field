using LinearAlgebra
using DelimitedFiles
using Statistics
#using Plots
using DelimitedFiles
include("Geometry.jl")

function get_H_pi_flux(R,inter_vector,N,t)  # attention: H_MF is on basis of sites, H_pi_flux is spin-degen
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
                if nearest(R[i],R2[j])
                    x=R[i][1]-R2[j][1]
                    y=R[i][2]-R2[j][2]
                    hop=0
                    if x*y<-0.01
                        hop=t
                    elseif x>0.1&&y>0.1
                        if mod(round(-R[i][2]*2/sqrt(3)),2)==0  # from B to A
                            hop=-t
                        else
                            hop=t
                        end
                    elseif x<-0.1&&y<-0.1
                        if mod(round(-R[i][2]*2/sqrt(3)),2)==0
                            hop=t
                        else
                            hop=-t
                        end
                    else
                        if mod(round(-R[i][2]*2/sqrt(3)),2)==0
                            hop=t
                        else
                            hop=-t
                        end
                    end
                    H[2*i-1,2*j-1]=H[2*i-1,2*j-1]+hop
                    H[2*i,2*j]=H[2*i,2*j]+hop
                end
            end
        end
    end
    return H
end

function get_H_uniform_hop(R,inter_vector,N,t)
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
                if nearest(R[i],R2[j])
                    H[2*i-1,2*j-1]=H[2*i-1,2*j-1]+t
                    H[2*i,2*j]=H[2*i,2*j]+t
                end
            end
        end
    end
    return H
end

function get_Chi_pi_flux(R,inter_vector,N,t)  
    Chi=zeros(Complex,N,N)
    n_inter=size(inter_vector,1)

    for nn=1:n_inter
        r=[]
        for ii=1:N
            push!(r,inter_vector[nn])
        end
        R2=R+r

        for i=1:N
            for j=1:N
                if nearest(R[i],R2[j])
                    x=R[i][1]-R2[j][1]
                    y=R[i][2]-R2[j][2]
                    hop=0
                    if x*y<-0.01
                        hop=t
                    elseif x>0.1&&y>0.1
                        if mod(round(-R[i][2]*2/sqrt(3)),2)==0  # from B to A
                            hop=-t
                        else
                            hop=t
                        end
                    elseif x<-0.1&&y<-0.1
                        if mod(round(-R[i][2]*2/sqrt(3)),2)==0
                            hop=t
                        else
                            hop=-t
                        end
                    else
                        if mod(round(-R[i][2]*2/sqrt(3)),2)==0
                            hop=t
                        else
                            hop=-t
                        end
                    end
                    Chi[i,j]=Chi[i,j]-hop        # hop in H has opposite sign to Chi
                end
            end
        end
    end
    return Chi
end

function get_Chi_uniform_hop(R,inter_vector,N,t)
    Chi=zeros(Complex,N,N)
    n_inter=size(inter_vector,1)

    for nn=1:n_inter
        r=[]
        for ii=1:N
            push!(r,inter_vector[nn])
        end
        R2=R+r

        for i=1:N
            for j=1:N
                if nearest(R[i],R2[j])
                    Chi[i,j]=Chi[i,j]-t      # hop in H has opposite sign to Chi
                end
            end
        end
    end
    return Chi
end
