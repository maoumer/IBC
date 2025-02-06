# Computes barrier certificate (standard vs IBC) using SOS for deterministic systems

# include important libraries
using JuMP
using MosekTools
using DynamicPolynomials
using MultivariatePolynomials
using LinearAlgebra
using TSSOS # important for SOS, see https://github.com/wangjie212/TSSOS

ϵ = 10^(-5) # B(x) ≥ ϵ instead of B(x) > 0
sos_tol = 1 # the maximum degree of unknown SOS polynomials = deg + sos_tol 
error = 5   # precision digit places

# polynomial barrier certificate of degree deg
function bc_standard(deg)
    # synthesize BC by using the standard formulation
    # deg: degree of BC template
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
    
    B, Bc, Bb = add_poly!(model, vars, deg) # generate polynomial template with given variables and degree
    Bf = B(vars=>f) #B(f(x)) https://juliapackages.com/p/dynamicpolynomials

    # see TSSOS https://github.com/wangjie212/TSSOS
    # 1st condition, initial set >= 0 by default (B ≤ 0i.e. -B ≥ 0)
    model,info1 = add_psatz!(model, -B, vars, go, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true) #-B >= 0
    # 2nd condition, unsafe set
    model,info2 = add_psatz!(model, B-ϵ , vars, gu, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    # 3rd and last condition, non-increasing
    model,info3 = add_psatz!(model, B-Bf, vars, g, [], div(maxdegree(B-Bf)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    
    optimize!(model) # solve for coefficients
    status = termination_status(model)
    Bc = value.(Bc)  # get the values of each coefficient
    for i = eachindex(Bc)
        Bc[i] = round(Bc[i]; digits = error) # round to order of error
    end
    # status might be optimal but if all Bc approx 10^{-error}, it's essentially 0 so OPTIMAL != BC.
    return (status,Bc'*Bb) # optimization status and barrier certificate function
end

# k = 1 -> bc_standard
function ibc(deg,k)
    # search IBC by using the given conditions
    # deg: degree of IBC template
    # k: max # BCs
    model = Model(optimizer_with_attributes(Mosek.Optimizer))
    set_optimizer_attribute(model, MOI.Silent(), true)
    
    # julia 1-indexing, i = 1 (i = 0 in paper)
    B, Bc, Bb = add_poly!(model, vars, deg) # generate polynomial template with given variables and degree
    Bf = B(vars=>f) #B(f(x)) https://juliapackages.com/p/dynamicpolynomials
    
    # 1st exclusive condition for B_0, >= 0 by default
    model,info1 = add_psatz!(model, -B, vars, go, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    # 2nd condition for B_0, unsafe states
    model,info2 = add_psatz!(model, B-ϵ , vars, gu, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    Bs  = [B] # store all unknown polynomial B(x) templates
    Bcs = [Bc] # store all corresponding unknown polynomial coefficients
    Bfs = [Bf] # store all corresponding unknown B(f(x)) templates
    for i = 2:k+1 # julia 1-indexing -> 2 <= i <= k+1 (1 <= i <= k in paper)
        B, Bc, Bb = add_poly!(model, vars, deg) 
        Bf = B(vars=>f)
        push!(Bs, B) # append to Bs in python
        push!(Bcs, Bc)
        push!(Bfs, Bf)
        # 2nd condition for B_i
        model,info2 = add_psatz!(model, B-ϵ , vars, gu, [], div(deg+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
        # 3rd condition for B_i(x) = Bs[i-1], B_(i+1)(f(x)) = Bf
        model,info3 = add_psatz!(model, Bs[i-1]-Bf, vars, g, [], div(maxdegree(Bs[i-1]-Bf)+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    end
    # 4th and last exclusive condition for B_k
    model,info3 = add_psatz!(model, Bs[k+1]-Bfs[k+1], vars, g, [], div(maxdegree(Bs[k+1]-Bfs[k+1])+sos_tol,2), QUIET=true, CS=false, TS=false, Groebnerbasis=true)
    
    optimize!(model) #solve for coefficients
    status = termination_status(model)  #all_variables(model)
    IBC = [] # store B_i
    # println(MosekTools.getprimalobj(model), MosekTools.getdualobj(model))
    for i = eachindex(Bcs)
        Bc = value.(Bcs[i]) # get the values of each coefficient for B_i
        for j = eachindex(Bc)
            Bc[j] = round(Bc[j]; digits = error) # round to order of error
        end
        push!(IBC, Bc'*Bb)
    end
    # status might be optimal but if all Bc approx 10^{-error}, it's essentially 0.
    return (status,IBC)
end


# Simulation
names = ["simple", "lotka-volt"]
system_data = ["vars: ","f: ","go: ","gu: ","g: "] # from systems files
max_deg, k_max = 5, 3
for name in names
    include("./systems/"*name*".jl"); # load system dynamics info
    file = open("./systems/"*name*"_system.txt", "w"); # open file to write/save system dynamics info
    for (j,k) in zip(system_data, [vars,f,go,gu,g]) # save these data in readable form
        write(file, j*"{")
        for i = 1:length(k)-1
            write(file, string(k[i])*", ")
        end
        write(file, string(last(k))*"}\n")
    end
    close(file)

    # print sufficient condition results for standard barrier certificate
    file = open("./systems/"*name*"_bc.txt", "w");
    for deg = 1:max_deg
        stats = @timed data = bc_standard(deg) # time execution
        status, B = data # get status and function
        write(file, "poly deg: "*string(deg)*"\n")
        write(file, "status: "*string(status)*"\n")
        write(file, Base.replace(string(B),"e"=>"*10^")*"\n")
        write(file, "time: "*string(stats.time)*"\n\n") 
    end
    close(file)

    # print IBC condition results
    file = open("./systems/"*name*"_ibc.txt", "w");
    for deg = 1:max_deg
        for k = 0:k_max
            stats = @timed data = ibc(deg, k)
            status, IBC = data
            write(file, "poly deg: "*string(deg)*",\t k: "*string(k)*"\n")
            write(file, "status: "*string(status)*"\n")
            write(file, Base.replace(string(IBC),"e"=>"*10^")*"\n")
            write(file, "time: "*string(stats.time)*"\n\n")
        end
    end
    close(file)
    println("Finished "*name)
end

