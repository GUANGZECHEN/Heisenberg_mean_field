using LinearAlgebra
using DelimitedFiles
using Statistics
#using Plots
using DelimitedFiles
include("Geometry.jl")
include("pi-flux.jl")
include("rearrange_results.jl")
include("matrix_initialization.jl")
include("get_H_MF.jl")

function get_new_parameters(Eigvals,Eigvecs,Fermi_level,N,R)
    M_new=zeros(Complex,2*N,2*N)
    Chi_new=zeros(Complex,N,N)
    n_new=zeros(Float64,N) # n_new is the filling on each site, which is contrained to 1

    n_bands=size(Eigvals,1)

    E_GS=0

    for ii=1:n_bands
        psi=Eigvecs[ii]
        if Eigvals[ii]+0.00001<Fermi_level
            filling=1
            E_GS=E_GS+Eigvals[ii]
        elseif Eigvals[ii]-0.00001>Fermi_level
            filling=0
        else
            filling=1/2
            E_GS=E_GS+1/2*Eigvals[ii]
        end
        
        for i=1:N
            M_new[2*i-1,2*i-1]=M_new[2*i-1,2*i-1]+(abs2(psi[2*i-1]))*filling
            M_new[2*i,2*i]=M_new[2*i,2*i]+(abs2(psi[2*i]))*filling
            M_new[2*i-1,2*i]=M_new[2*i-1,2*i]+(conj(psi[2*i-1])*psi[2*i])*filling
            M_new[2*i,2*i-1]=M_new[2*i,2*i-1]+(conj(psi[2*i])*psi[2*i-1])*filling

            n_new[i]=n_new[i]+(abs2(psi[2*i-1])+abs2(psi[2*i]))*filling

            
            for j=1:N
                #if nearest(R[i],R[j])   # hopping from nearest J1 interaction
                if j!=i
                    Chi_new[i,j]=Chi_new[i,j]+(conj(psi[2*i-1])*psi[2*j-1]+conj(psi[2*i])*psi[2*j])*filling   # in this way Chi is already complete and Hermitian
                end              
                #elseif next_nearest(R[i],R[j])   # hopping from next_nearest J2 interaction
                    #Chi_new[i,j]=Chi_new[i,j]+(conj(psi[2*i-1])*psi[2*j-1]+conj(psi[2*i])*psi[2*j])*filling
                #end
            end
            
        end

    end
    
    return M_new, Chi_new, n_new, E_GS

end

function MF_double_counting_energy(J1,J2,M,Chi,R,N)
    E=0
    for i=1:N
        for j=1:N
            if nearest(R[i],R[j])   # hopping from nearest J1 interaction
                E=E-(M[2*i-1,2*i-1]*M[2*j-1,2*j-1]+M[2*i-1,2*i]*M[2*j,2*j-1]+M[2*i,2*i-1]*M[2*j-1,2*j]+M[2*i,2*i]*M[2*j,2*j])*J1
                #println(E)
                E=E+Chi[i,j]*conj(Chi[i,j])*J1
                #println(E)              
            elseif next_nearest(R[i],R[j])   # hopping from next_nearest J2 interaction
                E=E-(M[2*i-1,2*i-1]*M[2*j-1,2*j-1]+M[2*i-1,2*i]*M[2*j,2*j-1]+M[2*i,2*i-1]*M[2*j-1,2*j]+M[2*i,2*i]*M[2*j,2*j])*J2
                E=E+Chi[i,j]*conj(Chi[i,j])*J2 
            end
        end
    end
    return E
end

function C3_Chi(Chi,N,R)
    Chi_new=zeros(Complex,N,N)
    for i=1:N
        for j=i:N
            Chi_new[i,j]=Chi[i,j]
            Chi_new[j,i]=conj(Chi[i,j])
        end
    end
    for i=1:N
        for k=i:N
            if abs(R[k][1]-((-1/2)*R[i][1]+(-sqrt(3)/2)*R[i][2]))+abs(R[k][2]-((sqrt(3)/2)*R[i][1]+(-1/2)*R[i][2]))<0.1
                for j=i+1:N
                    for l=k+1:N
                        if abs(R[l][1]-((-1/2)*R[j][1]+(-sqrt(3)/2)*R[j][2]))+abs(R[l][2]-((sqrt(3)/2)*R[j][1]+(-1/2)*R[j][2]))<0.1
                            Chi_new[k,l]=Chi[i,j]
                            Chi_new[l,k]=conj(Chi[i,j])
                        end
                    end                    
                end
            elseif abs(R[k][1]-((-1/2)*R[i][1]+(sqrt(3)/2)*R[i][2]))+abs(R[k][2]-((-sqrt(3)/2)*R[i][1]+(-1/2)*R[i][2]))<0.1
                for j=i+1:N
                    for l=k+1:N
                        if abs(R[l][1]-((-1/2)*R[j][1]+(sqrt(3)/2)*R[j][2]))+abs(R[l][2]-((-sqrt(3)/2)*R[j][1]+(-1/2)*R[j][2]))<0.1
                            Chi_new[k,l]=Chi[i,j]
                            Chi_new[l,k]=conj(Chi[i,j])
                        end
                    end                    
                end
            end
        end
    end
    return Chi_new
end

function add_magnetization_around_vacancy(M,m,vacancy_site)
    M[2*vacancy_site,2*vacancy_site]=0
    M[2*vacancy_site,2*vacancy_site-1]=im*10+1
    M[2*vacancy_site-1,2*vacancy_site]=-im*10+1
    M[2*vacancy_site-1,2*vacancy_site-1]=0
    M[2*vacancy_site-2,2*vacancy_site-2]=1
    M[2*vacancy_site-3,2*vacancy_site-3]=0
    M[2*(vacancy_site-m),2*(vacancy_site-m)]=0
    M[2*(vacancy_site-m)-1,2*(vacancy_site-m)-1]=1
    M[2*(vacancy_site-m+1),2*(vacancy_site-m+1)]=0
    M[2*(vacancy_site-m+1)-1,2*(vacancy_site-m+1)-1]=1
    M[2*(vacancy_site+m-2),2*(vacancy_site+m-2)]=1
    M[2*(vacancy_site+m-2)-1,2*(vacancy_site+m-2)-1]=0
    M[2*(vacancy_site+m-1),2*(vacancy_site+m-1)]=0
    M[2*(vacancy_site+m-1)-1,2*(vacancy_site+m-1)-1]=1
    return M
end

function solve_MF(Diff,acc,mixing,J1,J2,J3,N,M,Chi,Mu,R,inter_vector,lambda,lambda_2)
    Diff_M=zeros(Complex,2*N,2*N)
    Diff_Chi=zeros(Complex,N,N)
    Diff_filling=zeros(Float64,N)
    Eigvals=zeros(Float64,2*N)
    E_GS=0
    Fermi_level=0
    while Diff>acc
        M=M+mixing*Diff_M
        Chi=Chi+mixing*Diff_Chi
        
        #Chi=C3_Chi(Chi,N,R)
       
        for i=1:N
            Mu[i]=Mu[i]-0.1*mixing*tanh(Diff_filling[i])
        end

        Eigvals, Eigvecs=get_Syst(J1,J2,J3,N,M,Chi,Mu,R,inter_vector,lambda,lambda_2)
        Fermi_level=median(Eigvals)
        M_new, Chi_new, n_new, E_GS=get_new_parameters(Eigvals,Eigvecs,Fermi_level,N,R)
    
        Diff_M=M_new-M
        Diff_Chi=Chi_new-Chi
        Diff_filling=ones(Float64,N)-n_new
        Diff=norm(Diff_M)+norm(Diff_Chi)+norm(Diff_filling)

        println(Diff)

    end
    return Diff,Chi,Mu,M,Eigvals,E_GS,Fermi_level
end

function main()
    n=4
    m=4

    R,inter_vector=get_lattice("PBC","triangular",n,m)
    #R,inter_vector=get_lattice_honeycomb_nanoribbon(n,m)
    #R,inter_vector=get_lattice_triangular_nanoribbon(n,m)
    vacancy_site=6
    R=add_vacancy(R,vacancy_site)
    N=size(R,1)
    #print(R,N)

    filling=1/2   

    J1=1
    J2=0
    J3=0
    lambda=0.01
    lambda_2=0.99


    M=initialize_magnetization_matrix(N,"ferro_z")

    Chi=initialize_hop_matrix(R,inter_vector,"random",1)


    Mu=zeros(Complex,N)
    Eigvals=zeros(Float64,2*N)
    #Chi=C3_Chi(Chi,N,R)

    Diff=10

    Diff,Chi,Mu,M,Eigvals,E_GS,Fermi_level=solve_MF(Diff,1,0.1,J1,J2,J3,N,M,Chi,Mu,R,inter_vector,lambda,lambda_2)
    #Diff,Chi,Mu,M,Eigvals,E_GS,Fermi_level=solve_MF(Diff,0.5,0.001,J1,J2,J3,N,M,Chi,Mu,R,inter_vector,lambda,lambda_2)
    Diff,Chi,Mu,M,Eigvals,E_GS,Fermi_level=solve_MF(Diff,0.01,0.1,J1,J2,J3,N,M,Chi,Mu,R,inter_vector,lambda,lambda_2)

    net_mag_z=zeros(Float64,N)
    net_mag_y=zeros(Float64,N)
    net_mag_x=zeros(Float64,N)
    net_mag=zeros(Float64,N)
    Mu_ave=0
    for i=1:N
        net_mag_z[i]=abs(M[2*i-1,2*i-1])-abs(M[2*i,2*i])
        net_mag_y[i]=im*M[2*i,2*i-1]-im*M[2*i-1,2*i]
        net_mag_x[i]=M[2*i,2*i-1]+M[2*i-1,2*i]
        net_mag[i]=sqrt(net_mag_x[i]^2+net_mag_y[i]^2+net_mag_z[i]^2)
        Mu_ave=Mu_ave+Mu[i]/N
    end
    println("net_mag: ", net_mag)

    println(norm(Chi))
    #println("MF_double_counting_energy: ", MF_double_counting_energy(J1,J2,M,Chi,R,N))
    println("E_GS: ", E_GS-Mu_ave)
    #println("Fermi_level: ", Fermi_level)
    println("Chi: ",Chi)
    #println("M: ",M)
    println("Mu: ", Mu)
    
    #H_final=-1*get_H_MF(J1,J2,N,M,Chi,Mu,R,inter_vector)+2*J1*get_H_pi_flux_PBC(R,inter_vector,N,0.39)
    #M0=Matrix{Complex}(I,2*N,2*N)/2
    H_MF_final=lambda*get_H_MF(J1,J2,J3,N,M,Chi,Mu,R,inter_vector)+lambda_2*2*abs(J1)*get_H_pi_flux(R,inter_vector,N,1)+chemical_potential(N,Mu)
    println(H_MF_final)
    #println(inter_vector)

    results_n, results_nn=rearrange_results_hop(R,H_MF_final,inter_vector)
    results_mag=rearrange_results_mag(R,net_mag)
    results_mag_3D=rearrange_results_mag_3D(R,net_mag_x,net_mag_y,net_mag_z)
    writedlm( "size2_n_hop_J2=1_pyplot.csv",  results_n, ',')
    writedlm( "size2_nn_hop_J2=1_pyplot.csv",  results_nn, ',')
    writedlm( "sites_size2.csv",  results_mag, ',')
    writedlm( "reconstructed_spectra.csv", Eigvals, ',')
    writedlm( "mag_3D.csv",  results_mag_3D, ',')

end

function test()
    n=4
    m=2

    #R,inter_vector=get_lattice("OBC","triangular",n,m)
    R,inter_vector=get_lattice_honeycomb_nanoribbon(n,m)
    #R,inter_vector=get_lattice_triangular_nanoribbon(n,m)
    vacancy_site=6
    #R=add_vacancy(R,vacancy_site)
    N=size(R,1)
    #print(R,N)

    filling=1/2   

    J1=1
    J2=0
    J3=0


    M=initialize_magnetization_matrix(N,"random")

    Chi=initialize_hop_matrix(R,inter_vector,"uniform",1) 
    println("M: ",M)
    println("Chi: ",Chi)
    Mu=0

    H_MF_final=get_H_MF(J1,J2,J3,N,M,Chi,Mu,R,inter_vector)
    println(H_MF_final)
end



main()
#test()




