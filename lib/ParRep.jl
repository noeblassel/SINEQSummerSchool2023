module ParRep

    export GenParRepAlgorithm

    using Random, Base.Threads, Logging

    Base.@kwdef mutable struct GenParRepAlgorithm{S,P,M,K,R,X,L}
        N::Int #number of replicas

        # algorithm parameters
        simulator::S # a method to evolve the microscopic dynamics

        dephasing_checker::P # an object to check if replicas have dephased
        macrostate_checker::M # an object to check the macrostate
        replica_killer::K # an object to kill the replicas

        logger::L # log characteristics of the trajectory / debug 

        # Internal variables

        rng::R = Random.default_rng()
        reference_walker::X
        replicas::Vector{X} = typeof(reference_walker)[]

        n_initialisation_ticks::Int=0 # number of calls to simulator in the initialisation step
        n_dephasing_ticks::Int=0 # number of (parallel) calls to the simulator in the dephasing/decorrelation step
        n_parallel_ticks::Int=0 # number of (parallel) simulation steps in parallel step
        n_transitions::Int=0 # number of observed transitions
        simulation_time::Int=0 # simulation time -- in number of equivalent steps of the DNS Markov chain
        wallclock_time::Int=0
    end

    function simulate!(alg::GenParRepAlgorithm, n_transitions)
        current_macrostate = get_macrostate!(alg.macrostate_checker,alg.reference_walker,nothing)
        killed_ixs = Int[]

        while alg.n_transitions < n_transitions
            ## === INITIALISATION PHASE ===
            initialization_step = 0

            while current_macrostate === nothing
                # println("initialising --- iteration $(initialization_step)")
                update_microstate!(alg.reference_walker,alg.simulator)
                # println("\treference walker: ",alg.reference_walker)
                alg.n_simulation_ticks += 1
                current_macrostate = get_macrostate!(alg.macrostate_checker,alg.reference_walker,current_macrostate)
                initialization_step +=1
                log_state!(alg.logger,:initialization; algorithm = alg)
            end
            @info "Initialised in state $(current_macrostate)"

            alg.n_initialisation_ticks += initialization_step
            alg.simulation_time += initialization_step
            alg.wallclock_time += initialization_step

            ## === DECORRELATION/DEPHASING === 
            empty!(alg.replicas)
            for i=1:alg.N
                push!(alg.replicas,copy(alg.reference_walker))
            end
            
            has_dephased = false
            dephasing_step = 0

            while !has_dephased
                dephasing_step += 1
                # println("Dephasing --- iteration $(dephasing_step)")
                update_microstate!(alg.reference_walker,alg.simulator)
                ref_macrostate = get_macrostate!(alg.macrostate_checker,alg.reference_walker,current_macrostate)
                if check_death(alg.replica_killer,ref_macrostate,current_macrostate,alg.rng)
                log_state!(alg.logger,:transition; algorithm = alg,current_macrostate=current_macrostate,new_macrostate=ref_macrostate,exit_time=0,dephasing_time=dephasing_step)
                     current_macrostate = ref_macrostate
                    # println("\tReference walker has crossed to state $(current_macrostate)")
                    # print("Q")
                    alg.n_transitions += 1
                    break
                end
                # println("\treference walker: ",alg.reference_walker)
                empty!(killed_ixs)
                # check if reference_walker has escaped


                @threads for i=1:alg.N # to parallelize in production
                    update_microstate!(alg.replicas[i],alg.simulator)
                    rep_macrostate = get_macrostate!(alg.macrostate_checker,alg.replicas[i],current_macrostate)

                    if check_death(alg.replica_killer,rep_macrostate,current_macrostate,alg.rng)
                        push!(killed_ixs,i)
                        # println("\t\tReplica $i crossed to state $(rep_macrostate)")
                    end

                end
                log_state!(alg.logger,:dephasing; algorithm = alg,killed_ixs=killed_ixs)
                survivors = setdiff(1:alg.N,killed_ixs)

                # println("\t$(length(killed_ixs)) replicas have been killed")

                if length(survivors)==0
                    throw(DomainError("Extinction event!"))
                end
                
                @threads for i=killed_ixs
                    alg.replicas[i] = copy(alg.replicas[rand(alg.rng,survivors)])
                end
                
                has_dephased = check_dephasing!(alg.dephasing_checker,alg.replicas,current_macrostate,dephasing_step)

                # Obar = sum(alg.dephasing_checker.means;dims = 2) / alg.N
                # numerator = sum(@. (alg.dephasing_checker.sq_means -2alg.dephasing_checker.means*Obar + Obar^2);dims=2)
                # denominator = sum(alg.dephasing_checker.sq_means - alg.dephasing_checker.means .^ 2;dims=2)
                # println("\tGelman-Rubin diagnoses: ",numerator./denominator)
            end

            alg.n_dephasing_ticks += dephasing_step
            alg.simulation_time += dephasing_step
            alg.wallclock_time += dephasing_step

            if has_dephased
                @info "Successfully dephased in state $current_macrostate"
                ## === PARALLEL PHASE === 
                killed = false
                i_min = nothing
                new_macrostate = nothing
                parallel_step = 0
                
                while !killed
                    @threads for i=1:alg.N
                        if killed == true
                            break
                        end

                        update_microstate!(alg.replicas[i],alg.simulator)
                        rep_macrostate = get_macrostate!(alg.macrostate_checker,alg.replicas[i],current_macrostate)
                        
                        if check_death(alg.replica_killer,rep_macrostate,current_macrostate,alg.rng)
                            killed = true
                            i_min = i
                            new_macrostate = rep_macrostate
                            @info "Replica $i escaped to state $(rep_macrostate)"
                            break
                        end
                    end

                    parallel_step += 1
                    log_state!(alg.logger,:parallel; algorithm = alg)
                end

                
                alg.reference_walker = copy(alg.replicas[i_min])
                exit_time = (alg.N*parallel_step + i_min)

                log_state!(alg.logger,:transition; algorithm = alg,current_macrostate=current_macrostate,new_macrostate=new_macrostate,exit_time=exit_time,dephasing_time=dephasing_step)
                
                alg.n_parallel_ticks += parallel_step
                alg.wallclock_time += parallel_step
                alg.n_transitions += 1
                alg.simulation_time += exit_time
                current_macrostate = new_macrostate
            end
        end
    end


    # dummy method definitions

    function check_dephasing!(checker,replicas,current_macrostate,step_n) end
    function get_macrostate!(checker,walker,current_macrostate) end
    function update_microstate!(simulator,walker) end
    function check_death(checker,macrostate_a,macrostate_b,rng) end
    function log_state!(logger,event; kwargs...) end
end