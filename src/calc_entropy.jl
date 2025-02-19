#!/usr/bin/env julia
# include("Functions.jl")
# include("lammpstrj_parser.jl")
# include("xyz_parser.jl")
# include("main.jl")

using EntMix

using ArgParse
using Base.Threads
using DelimitedFiles
using Pkg

argparser = ArgParseSettings()
@add_arg_table! argparser begin
    "--radial_factor", "-f"
        help="Factor to multiply the VdW radii of the atoms"
        default=.6
        arg_type=Float64
    "--maxstep", "-m"
        help="Maximum step to use (default all steps)"
        default=nothing
        arg_type=Int
    "--startstep", "-s"
        help="Initial step to use default to first step in file"
        default=nothing
        arg_type=Int
    "--nthstep" 
        help="Only use every nth step"
        default=nothing
        arg_type=Int
    "--densfunc", "-d" 
        help="Function to use for smearing. Available options are: slater, gaus"
        default="slater"
        arg_type=String
    "--nonperiodic", "-n"
        help="If set, the box will be treated as non periodic. Only applicable for lammpstrj files since xyz files are always non periodic."
        action="store_true"
    "--molats1"
        help="Number of atoms in the first molecule type. If set the program will normalize the density by this amount."
        arg_type=Int
        default=1
    "--molats2"
        help="Number of atoms in the second molecule type. If set the program will normalize the density by this amount."
        arg_type=Int
        default=1
    "--outfile", "-o"
        help="Output file to save the entropy values, if not set, the values will be printed to the console.
Only applicable for lammpstrj files. If set threads can be used with julia -t <nthreads>."
        arg_type=String
    "--debug", "-v"
        help="Print debug information"
        action="store_true"
    "trajfile"
        help="xyz or lammpstrj file to use. 
Using a xyz file will result in a non periodic box calculation. 
Lammpstrj files will be treated as periodic if not defined otherwise with --nonperiodic."
        required=true
        arg_type=String
    "amatoms"
        help="Total number of type 1 atoms in the system. They must be the first atoms in the file."
        required=true
        arg_type=Int
end

args = parse_args(argparser)
if args["debug"]
    ENV["JULIA_DEBUG"] = "EntMix,Main"
end
radial_factor = args["radial_factor"]
file = args["trajfile"]
ftype = split(file, ".")[end]
if ftype == "lammpstrj" ||  ftype == "lmp"
    trajtype = "lammpstrj"
elseif ftype == "xyz"
    trajtype = "xyz"
else
    println("Invalid file type")
    println("Only .xyz and .lammpstrj files are supported")
    exit(1)
end
n_atoms = args["amatoms"]
maxstep = args["maxstep"]
startstep = args["startstep"]
outfile = args["outfile"]
nth = args["nthstep"]
if args["densfunc"] == "slater"
    dfunc = EntMix.slater
elseif args["densfunc"] == "gaus"
    dfunc = EntMix.gaus
# elseif args["densfunc"] == "rect"
#     dfunc = EntMix.rect
else
    println("Invalid density function")
    exit(1)
end
molnorm = args["molats1"], args["molats2"]
@debug "Molecule atoms: $molnorm"

pinfo = Pkg.dependencies()[Pkg.project().dependencies["EntMix"]]
@info "Running EntMix.jl version: $(pinfo.version)"
if pinfo.git_source == nothing
    # if added in dev mode
    moduledir = pinfo.source
else
    # if added via Pkg.add
    moduledir = pinfo.git_source
end

moduledir 
gitstatus = strip(read(`git -C $moduledir status --porcelain`, String))
@info "Running at git commit: $(read(`git -C $moduledir rev-parse HEAD`, String))"
if gitstatus != ""
    @warn "The module directory is not clean. Please commit your changes."
    println("Current Diff:", read(`git -C $moduledir diff`, String))
end

# run Main function
if trajtype == "lammpstrj"
    EntMix.lammpstrj_entropy(
                             file, 
                             n_atoms; 
                             maxstep=maxstep, 
                             startstep=startstep, 
                             radial_factor=radial_factor, 
                             dfunc=dfunc, 
                             outfile=outfile, 
                             na1=molnorm[1], na2=molnorm[2], 
                             nth=nth)
elseif trajtype == "xyz"
    EntMix.xyz_entropy(file, n_atoms, maxstep, startstep; radial_factor=radial_factor, dfunc=dfunc, na1=molnorm[1], na2=molnorm[2])
end

