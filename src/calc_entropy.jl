#!/usr/bin/env julia

using EntMix 
import Chemfiles
using ArgParse
using Base.Threads
# using DelimitedFiles

argparser = ArgParseSettings()
@add_arg_table! argparser begin
    "--debug", "-v"
        help="Print debug information"
        action="store_true"
    "--startstep"
        help="Initial step to use default to first step in file"
        default=0
        arg_type=Int
    "--endstep"
        help="Final step to use (default last step in file)"
        default=999999999999999
        arg_type=Int
    "--stepinterval"
        help="Interval between steps to use (default 1, using each step in the file)"
        default=1
        arg_type=Int
    "--sigma"
        help="Width of the Gaussian kernel used for entropy calculation (default 1.0) unit of length in trajectory"
        default=1.0
        arg_type=Float64
    "--sigmatype"
        help="Type of sigma to use: *homo, VDW, Covalent"
        default="homo"
        arg_type=String
    "--outfile", "-o"
        help="Output file to write results to (default stdout)"
        default=""
        arg_type=String
    "--append", "-a"
        help="Append to outfile"
        action="store_true"
    "trajfile"
        help="trajectory file to use."
        required=true
        arg_type=String
end

args = parse_args(argparser)
if args["debug"]
    ENV["JULIA_DEBUG"] = "EntMix,Main"
end

function main()
    trajectory = Chemfiles.Trajectory(args["trajfile"])
    firstframe = read(trajectory)
    if Chemfiles.type(firstframe[1]) == "" 
        Molecule.type_from_name!(firstframe)
    end
    mols = Molecule.get_molecules(firstframe)
    moltypes = Molecule.mol_types(firstframe, mols)
    @info """Detected $(length(moltypes)) molecule types:
        $(join(keys(moltypes), "\n"))"""

    @info "Calculating mixing entropy between these molecules"
    maxstep = length(trajectory) -1
    if args["endstep"] < maxstep 
        maxstep = args["endstep"] 
    end
    if args["startstep"] > 0 
        startstep = args["startstep"]
    else
        startstep = 0
    end
    entropies = Vector{Vector{Float64}}()
    traj_lock = ReentrantLock()
    entro_lock = ReentrantLock()
    frame = nothing
    @threads for stepid in startstep:args["stepinterval"]:maxstep
        @lock traj_lock frame = Chemfiles.read_step(trajectory, stepid)
        if args["sigmatype"] != "homo"
            if Chemfiles.type(frame[1]) == ""
                Molecule.type_from_name!(frame)
            end
        end
        atom_coll = [reduce(vcat, [mols[mid] for mid in moltypes[key]]) for key in keys(moltypes)]
        entropy_val = [stepid, EntMix.entropy(
                                              frame, 
                                              atom_coll, 
                                              args["sigma"]; 
                                              baselength=args["sigmatype"]
                                             )]
        @lock entro_lock push!(entropies, entropy_val)
        @debug entropies[end]
    end
    if args["outfile"] != ""
        if args["append"]
            open(args["outfile"], "a") do io
            for ent in sort(entropies, by=x->x[1])
                println(io, join(ent, ", "))
            end
            end
        else
            open(args["outfile"], "w") do io
            for ent in sort(entropies, by=x->x[1])
                println(io, join(ent, ", "))
            end
            end
        end
        # writedlm(args["outfile"], sort(entropies, by=x->x[1]))
    else
        for ent in sort(entropies, by=x->x[1])
            println(join(ent, ", "))
        end
    end
end
@time main()

