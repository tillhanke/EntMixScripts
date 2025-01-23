using EntMix
using Base.Threads
using ArgParse

include("./cube.jl")

argparser = ArgParseSettings()
@add_arg_table! argparser begin 
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
    "--out", "-o"
        help="Outputfile (default [trajfile].cube)"
        required=false
        arg_type=String
        default=nothing
    "--timestep", "-t"
        help="Timestep, to be used (defaults to first in file)"
        required=false
        arg_type=Int
        default=0
     "--gridpoints", "-g"
        help="Number of gridpoints in each direction (default 60)"
        required=false
        arg_type=Int
        default= 60  
    "--radial_factor", "-f"
        help="Factor to multiply the VdW radii of the atoms"
        default=.6
        arg_type=Float64
    "-v"
        help="Verbose output"
        action=:store_true
end

args = parse_args(argparser)
if args["v"]
    ENV["JULIA_DEBUG"] = "Main,EntMix"
end
@debug args

if !(endswith(args["trajfile"], ".lammpstrj") || endswith(args["trajfile"], ".lmp"))
    error("Only lammpstrj files are supported")
end

if args["timestep"] == 0
    elems, coords, boxs = EntMix.parse_lammpstrj(
        args["trajfile"];
        start=args["timestep"],
        amsteps=1    
       )
else
    elems, coords, boxs = EntMix.parse_lammpstrj(
        args["trajfile"];
        start=args["timestep"],
        stop=args["timestep"]
    )
end

box = boxs[1, :, :]
elems = elems[1, :]
coords = coords[1, :, :]

# move box to origin
if box[1,:] != [0, 0, 0]
    box -= box[1,:]
    box[1,:] = [0, 0, 0]
end
cubebox = [
           box[2,1] 0 0
           0 box[2,2] 0
           0 0 box[2,3]
          ]
cubedefinition = (
                  box = cubebox./args["gridpoints"], 
                  repetitions = [args["gridpoints"], args["gridpoints"], args["gridpoints"]], 
                  origin = [0, 0, 0]
                 )

cubegrid = create_grid(cubedefinition...)
densfrac = Array{Float64, 3}(undef, args["gridpoints"], args["gridpoints"], args["gridpoints"])
EntMix.wrap!(coords, box)
cA, elA = EntMix.add_periodic!(coords[1:args["amatoms"], :], elems[1:args["amatoms"], :], box; delta=10)
@debug "amount added: $(size(cA))"
cA = [
    cA
    coords[1:args["amatoms"], :]
   ]
cB, elB = EntMix.add_periodic!(coords[args["amatoms"]+1:end, :], elems[args["amatoms"]+1:end, :], box; delta=10)
@debug "amount added: $(size(cB))"
cB = [
    cB
    coords[args["amatoms"] + 1:end, :]
   ]
# coords, elems = EntMix.add_periodic!(coords[args["amatoms"]+1:end, :], elems, box; delta=10)

sigA = [EntMix.VDWradii[elem]*args["radial_factor"] for elem in elA]
sigB = [EntMix.VDWradii[elem]*args["radial_factor"] for elem in elB]

for i in 1:args["gridpoints"]
    @debug "Calculating xi=$i"
    for j in 1:args["gridpoints"]
        @threads for k in 1:args["gridpoints"]
            densfrac[i, j, k] = EntMix.density_fraction(
                EntMix.dens(
                    cubegrid[i, j, k],
                    atoms=cA,
                    sigma=sigA
                ),
                EntMix.dens(
                    cubegrid[i, j, k],
                    atoms=cB,
                    sigma=sigB
                ),
            )
        end
    end
end

outfile = args["out"] == nothing ? replace(args["trajfile"], r"\.[^\.]*$" => ".cube") : args["out"]
@debug size(densfrac)

write_cube(outfile, [EntMix.atomic_numbers[elem] for elem in elems], coords, cubedefinition..., densfrac; inputunits="ANG")

