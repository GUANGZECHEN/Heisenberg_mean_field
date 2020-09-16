using LinearAlgebra
using DelimitedFiles
using Statistics
#using Plots
using DelimitedFiles

function get_lattice(BC,geometry,n,m)
    if geometry=="triangular"
        a1=[1,0]
        a2=[1/2,-sqrt(3)/2]
        R=geometry_simple(n,m,a1,a2)
    elseif geometry=="square"
        a1=[1,0]
        a2=[0,1]
        R=geometry_simple(n,m,a1,a2)
    elseif geometry=="honeycomb"
        a1=[sqrt(3),0]
        a2=[sqrt(3)/2,-3/2]
        a3=[sqrt(3)/2,-1/2]
        R=geometry_bipartite(n,m,a1,a2,a3)
    else
        println("invalid geometry")
    end

    if BC=="PBC"
        inter_vector=get_inter_vector(n,m,a1,a2)
    elseif BC=="OBC"
        inter_vector=[[0,0]]
    else
        println("invalid boundary condition")
    end
    
    return R,inter_vector
end

function geometry_simple(n,m,a1,a2)
    R=[]
    for i=0:n-1
        for j=0:m-1
            push!(R,i*a1+j*a2)
        end
    end
    return R
end

function geometry_bipartite(n,m,a1,a2,a3)
    R=[]
    for i=0:n-1
        for j=0:m-1
            for k=0:1
                push!(R,i*a1+j*a2+k*a3)
            end
        end
    end
    return R
end

function get_inter_vector(n,m,a1,a2)
    inter_vector=[[0,0],n*a1,-n*a1,m*a2,-m*a2,m*a2-n*a1,n*a1-m*a2,m*a2+n*a1,-n*a1-m*a2]
    return inter_vector
end

function get_lattice_honeycomb_nanoribbon(n,m)
    a1=[sqrt(3),0]
    a2=[0,-3]
    b1=[sqrt(3)/2,-1/2]
    b2=[sqrt(3)/2,-3/2]
    b3=[0,-2]

    R=[]
    for i=0:n-1
        for j=0:m-1
            push!(R,i*a1+j*a2)
            push!(R,i*a1+j*a2+b1)
            push!(R,i*a1+j*a2+b2)
            push!(R,i*a1+j*a2+b3)
        end
    end
    inter_vector=[[0,0],n*a1,-n*a1]
    return R,inter_vector
end

function get_lattice_triangular_nanoribbon(n,m)
    a1=[1,0]
    a2=[0,-sqrt(3)]
    b1=[1/2,-sqrt(3)/2]


    R=[]
    for i=0:n-1
        for j=0:m-1
            push!(R,i*a1+j*a2)
            push!(R,i*a1+j*a2+b1)
        end
    end
    inter_vector=[[0,0],m*a2,-m*a2]
    return R,inter_vector
end

function geometry_triangular(n,a1,a2)
    R=[]
    for i=0:n-1
        for j=0:n-1-i
            push!(R,i*a1+j*a2)
        end
    end
    return R
end

function geometry_zigzag(size)
    a1=[1,0]
    a2=[-1/2,sqrt(3)/2]
    a3=[-1/2,-sqrt(3)/2]
    R=[]
    for i=0:size
        for j=0:size
            if abs(i)+abs(j)<size+2
                push!(R,i*a1+j*a2)
            end
        end

        for l=1:size
            if abs(i)+abs(l)<size+2
                push!(R,i*a1+l*a3)
            end
        end
    end

    for j=1:size
        for l=1:size
            if abs(j)+abs(l)<size+2
                push!(R,l*a3+j*a2)
            end
        end
    end
    return R
end

function nearest(site_a,site_b)
    return 0.1<norm(site_a-site_b)<1.1
end

function next_nearest(site_a,site_b)
    return 1.1<norm(site_a-site_b)<1.9
end


function next_next_nearest(site_a,site_b)
    return 1.9<norm(site_a-site_b)<2.1
end

function check_R(R)
    N=size(R,1)
    for i=1:N
        for j=1:N
            if nearest(R[i],R[j])   
                println("near",R[i],R[j])
            elseif next_nearest(R[i],R[j])  
                println("n-near",R[i],R[j])
            end
        end
    end
end

function add_vacancy(R,vacant_site)
    deleteat!(R,vacant_site)
    return R
end





