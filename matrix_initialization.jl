using LinearAlgebra
#using DelimitedFiles
using Random
include("pi-flux.jl")

function initialize_magnetization_matrix(N,mode)
    if mode=="ferro_z"
        M=zeros(Complex,2*N,2*N)
        for i=1:N
            M[2*i-1,2*i-1]=1
        end
    elseif mode=="ferro_y"
        M=zeros(Complex,2*N,2*N)
        for i=1:N
            M[2*i-1,2*i]=im/2*(-1)
            M[2*i,2*i-1]=-im/2*(-1)
        end
    elseif mode=="ferro_x"
        M=zeros(Complex,2*N,2*N)
        for i=1:N
            M[2*i-1,2*i]=1/2
            M[2*i,2*i-1]=1/2
        end
    elseif mode=="non-magnetic"
        M=Matrix{Complex}(I,2*N,2*N)/2
    elseif mode=="random"
        M=random_Hermitian_matrix(2*N)
    else
        println("invalid magnetization mode")
    end
    return M
end

function initialize_hop_matrix(R,inter_vector,mode,t)  # R:geometry, inter_vector: PBC or OBC, t:magnitude of hopping
    N=size(R,1)
        
    if mode=="pi-flux"      
        Chi=get_Chi_pi_flux(R,inter_vector,N,t)
    elseif mode=="uniform"
        Chi=get_Chi_uniform_hop(R,inter_vector,N,t)
    elseif mode=="random"
        Chi=random_Hermitian_matrix(N)*t
    elseif mode=="no-hopping"
        Chi=zeros(Complex,N,N)
    else
        println("invalid hopping mode")
    end

    return Chi
end

function get_AF_triangular_nanoribbon(R,N)
    M=zeros(Complex,2*N,2*N)
    for i=1:N
        if mod(round(-R[i][2]*2/sqrt(3)),2)==0
            M[2*i-1,2*i-1]=1
        else
            M[2*i,2*i]=1
        end
    end
    return M
end

function lower_triangular_matrix(N) # generate lower triangular matrix of dimension N
    M=[i>j ? 1 : 0 for i=1:N, j=1:N]
    return M
end

function upper_triangular_matrix(N) # generate lower triangular matrix of dimension N
    M=[i<j ? 1 : 0 for i=1:N, j=1:N]
    return M
end

function random_Hermitian_matrix(N)
    M=rand(ComplexF64,(N,N))-ones(ComplexF64,(N,N))/2
    M=(M+adjoint(M))/2
    return M
end





