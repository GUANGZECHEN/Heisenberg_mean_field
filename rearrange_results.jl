using LinearAlgebra
using DelimitedFiles
include("Geometry.jl")

function rearrange_results(R,Chi)
    results_n=[]
    results_nn=[]
    N=size(R,1)
    for i=1:N
        for j=i:N
            if nearest(R[i],R[j])
                x=(R[i]+R[j])/2
                push!(x,norm(Chi[i,j]))
                push!(results_n,x)
            end
            if next_nearest(R[i],R[j])
                x=(R[i]+R[j])/2
                push!(x,norm(Chi[i,j]))
                push!(results_nn,x)
            end
        end
    end
    return results_n, results_nn    
end

function rearrange_results2(R,Chi)
    results_n=[]
    results_nn=[]
    N=size(R,1)
    for i=1:N
        for j=i:N
            if nearest(R[i],R[j])
                x=[]
                x=push!(x,R[i][1],R[i][2],R[j][1],R[j][2],norm(Chi[i,j]))
                push!(results_n,x)
            end
            if next_nearest(R[i],R[j])
                x=[]
                x=push!(x,R[i][1],R[i][2],R[j][1],R[j][2],norm(Chi[i,j]))
                push!(results_nn,x)
            end
        end
    end
    return results_n, results_nn    
end

function rearrange_results_hop(R,H,inter_vector)
    n_inter=size(inter_vector,1)
    results_positive=[]
    results_negative=[]
    N=size(R,1)

    for nn=1:n_inter
        #println(inter_vector[nn])
        r=[]
        for ii=1:N
            push!(r,inter_vector[nn])
        end
        R2=R+r

        for i=1:N
            for j=1:N
                hop=H[2*i,2*j]
                if nearest(R[i],R2[j]) || next_nearest(R[i],R2[j]) || next_next_nearest(R[i],R2[j])
                    #println(R[i],R2[j],hop)
                    x=[]
                    x=push!(x,R[i][1],R[i][2],R2[j][1],R2[j][2],abs(hop))
                    if real(hop)>0
                        push!(results_positive,x)
                    else
                        push!(results_negative,x)
                    end
                end
            end
        end
    end
    return results_positive, results_negative    
end

function rearrange_results_mag(R,mag)
    results=[]

    N=size(R,1)
    for i=1:N
        x=[]
        x=push!(x,R[i][1],R[i][2],mag[i])
        push!(results,x)
    end
    return results   
end

function rearrange_results_mag_3D(R,mag_x,mag_y,mag_z)
    results=[]

    N=size(R,1)
    for i=1:N
        x=[]
        x=push!(x,R[i][1],R[i][2],0,mag_x[i],mag_y[i],mag_z[i])
        push!(results,x)
    end
    return results   
end

