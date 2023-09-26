module ParRep

    export GenParRepAlgorithm

    using Random, Base.Threads, Logging

    Base.@kwdef mutable struct GenParRepAlgorithm{S,P,M,K,R,X,L}

        simulator::S
        dephasing_checker::P
        macrostate_checker::M
        replica_killer::K

        logger::L

        N::Int

        rng::R = Random.default_rng()
        reference_walker::X
        replicas::Vector{X} = typeof(reference_walker)[]

        n_initialisation_ticks::Int=0
        n_dephasing_ticks::Int=0
        n_parallel_ticks::Int=0

        n_transitions::Int=0
        simulation_time::Int=0
        wallclock_time::Int=0
    end

    function simulate!(alg::GenParRepAlgorithm, n_transitions)

        current_macrostate = get_macrostate!(alg.macrostate_checker,alg.reference_walker,nothing,0)
        alg.replicas = [copy(reference_walker) for i=1:alg.N]

        killed_ixs = Int[]

        while alg.n_transitions < n_transitions
            
            # Initialization step

            initialization_step = 0

            while current_macrostate === nothing
                initialization_step +=1
                update_microstate!(alg.reference_walker,alg.simulator)
                alg.n_simulation_ticks += 1
                current_macrostate = get_macrostate!(alg.macrostate_checker,alg.reference_walker,current_macrostate,initialization_step)
                log_state!(alg.logger,:initialization; algorithm = alg)
            end
            @info "Initialised in state $(current_macrostate)"

            alg.n_initialisation_ticks += initialization_step
            alg.simulation_time += initialization_step
            alg.wallclock_time += initialization_step

            # Decorrelation / dephasing step

            for i=1:alg.N
                push!(alg.replicas,copy(alg.reference_walker))
            end
            
            has_dephased = false
            dephasing_step = 0

            while !has_dephased
                dephasing_step += 1
                update_microstate!(alg.reference_walker,alg.simulator)
                ref_macrostate = get_macrostate!(alg.macrostate_checker,alg.reference_walker,current_macrostate,dephasing_step)

                if check_death(alg.replica_killer,ref_macrostate,current_macrostate,alg.rng)

                    # Failed decorrelation
                    log_state!(alg.logger,:transition; algorithm = alg,current_macrostate=current_macrostate,new_macrostate=ref_macrostate,exit_time=0,dephasing_time=dephasing_step)
                    current_macrostate = ref_macrostate
                    alg.n_transitions += 1
                    break
                end

                empty!(killed_ixs)

                @threads for i=1:alg.N
                    update_microstate!(alg.replicas[i],alg.simulator)
                    rep_macrostate = get_macrostate!(alg.macrostate_checker,alg.replicas[i],current_macrostate,dephasing_step)

                    if check_death(alg.replica_killer,rep_macrostate,current_macrostate,alg.rng)
                        push!(killed_ixs,i)
                    end

                end
                log_state!(alg.logger,:dephasing; algorithm = alg,killed_ixs=killed_ixs)
                survivors = setdiff(1:alg.N,killed_ixs)

                if length(survivors)==0
                    throw(DomainError("Extinction event!")) # All replicas have died
                end
                
                @threads for i=killed_ixs
                    alg.replicas[i] = copy(alg.replicas[rand(alg.rng,survivors)])
                end
                
                has_dephased = check_dephasing!(alg.dephasing_checker,alg.replicas,current_macrostate,dephasing_step)
            end

            alg.n_dephasing_ticks += dephasing_step
            alg.simulation_time += dephasing_step
            alg.wallclock_time += dephasing_step

            if has_dephased
                # Succesful decorrelation => Parallel exit phase
                @info "Successfully dephased in state $current_macrostate"

                killed = false
                i_min = nothing
                new_macrostate = nothing
                parallel_step = 0
                
                while !killed
                    parallel_step += 1
                    @threads for i=1:alg.N
                        if killed == true
                            break
                        end

                        update_microstate!(alg.replicas[i],alg.simulator)
                        rep_macrostate = get_macrostate!(alg.macrostate_checker,alg.replicas[i],current_macrostate,parallel_step)
                        
                        if check_death(alg.replica_killer,rep_macrostate,current_macrostate,alg.rng)
                            killed = true
                            i_min = i
                            new_macrostate = rep_macrostate
                            @info "Replica $i escaped to state $(rep_macrostate)"
                            break
                        end
                    end

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

        return alg.simulation_time,alg.wallclock_time
    end


    # fallback method definitions

    function check_dephasing!(checker,replicas,current_macrostate,step_n) end
    function get_macrostate!(checker,walker,current_macrostate,step_n) end
    function update_microstate!(simulator,walker) end
    function check_death(checker,macrostate_a,macrostate_b,rng) end
    function log_state!(logger,event; kwargs...) end
end