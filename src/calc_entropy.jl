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
    "--smearfunc"
        help="The type of smearing function to be used for density calculations. \nOptions are: gaus,slater,const,linear,epanechnikov,cubicspline,wendland,gengaus,cosine"
        default="gaus"
        arg_type=String
    "--outfile", "-o"
        help="Output file to write results to (default stdout)"
        default=""
        arg_type=String
    "--boxlengths", "-b"
        help="The lengths of your periodic box (e.g. 2.5,1.2,6.3). This code only uses orthogonal boxes. only used if not present in trajfile"
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

function parse_func(func_string)
    if func_string == "slater"
        return EntMix.slater
    elseif func_string == "gaus"
        return EntMix.gaus
    elseif func_string == "const"
        return EntMix.constant
    elseif func_string == "linear"
        return EntMix.linear
    elseif func_string == "epanechnikov"
        return EntMix.epanechnikov
    elseif func_string == "cubicspline"
        return EntMix.cubicspline
    elseif func_string == "wendland"
        return EntMix.wendland
    elseif func_string == "gengaus"
        return EntMix.gengauss
    elseif func_string == "cosine"
        return EntMix.cosine
    else
        @error func_string " is not a valid smearing function"
        exit(2)
    end
end

function main()
    trajectory = Chemfiles.Trajectory(args["trajfile"])
    firstframe = read(trajectory)
    c_box = args["boxlengths"]
    if c_box != "" 
        if Chemfiles.lengths(Chemfiles.UnitCell(firstframe)) ==  [0., 0., 0.]
            box = parse.(Float64, split(c_box, ','))
            cell = Chemfiles.UnitCell(box)
            Chemfiles.set_cell!(trajectory, cell)
        end
    end

    if Chemfiles.type(firstframe[1]) == "" 
        Molecule.type_from_name!(firstframe)
    end
    mols = Molecule.get_molecules(firstframe)
    moltypes = Molecule.mol_types(firstframe, mols)
    natoms = [length(mols[moltypes[key][1]]) for key in keys(moltypes)]
    @debug "Number of atoms per molecule type: $(natoms)"
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
    # traj_lock = ReentrantLock()
    entro_lock = ReentrantLock()
    # frame = nothing
    stepids = startstep:args["stepinterval"]:maxstep
    @debug "Reading the following stepids (chemfilenotation): $(stepids)"
    frames = [Chemfiles.read_step(trajectory, stepid) for stepid in stepids]
    @debug "type of first frames entry: $(typeof(frames[1]))"
    @debug "Box size of first frame: $(Chemfiles.lengths(Chemfiles.UnitCell(frames[1])))"
    if prod(Chemfiles.lengths(Chemfiles.UnitCell(frames[1]))) == 0
        @error "No Unit Cell is set for the first frame"
        exit(2)
    end

    @threads for id in eachindex(frames) 
        # @lock traj_lock frame = Chemfiles.read_step(trajectory, stepid)
        frame = frames[id]
        stepid = stepids[id]
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
                                              baselength=args["sigmatype"],
                                              dfunc=parse_func(args["smearfunc"]),
                                              natoms=natoms
                                             )]
        @debug "first atom position: $(Chemfiles.positions(frame)[:,1])"
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

