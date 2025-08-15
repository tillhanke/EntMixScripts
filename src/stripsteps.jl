#!/usr/bin/env julia
#
# this script can be used to strip certain steps from a trajectory file 
#
include("lammpstrj_parser.jl")
using ArgParse

function write_step(file, outfile)
    head = read_lammpstrj_head(file)
    global ffstep 
    global dts
    approx_b = position(trajfile)÷((head["timestep"]-ffstep)÷dts)
    for i in 1:head["n_atoms"]+9
        write(outfile, readline(file) * "\n")   
    end
    return approx_b
end

argparser = ArgParseSettings()
@add_arg_table! argparser begin
    "--laststep", "-l"
        help="last step to save defaults to the last step in the file"
        default=nothing
        arg_type=Int
    "--startstep", "-s"
        help="First step to save"
        default=nothing
        arg_type=Int
    "--skip", "-k"
        help="Number of steps to skip inbetween. (e.g. if 5 then only start, start+5, start+10, ...) will be written to output"
        default=nothing
        arg_type=Int    
    "--nsteps", "-n"
        help="Number of steps to write, starting from startstep. Overwrites --laststep. If not given, all steps until --laststep will be written"
        default=nothing
        arg_type=Int
    "--outfile", "-o"
        help="Output file, to save the stripped trajectory"
    "--debug", "-v"
        help="Print debug information"
        action="store_true"
    "trajfile"
        help="input trajectory file"
        required=true
        arg_type=String
end    

args = parse_args(argparser)
if args["debug"]
    ENV["JULIA_DEBUG"] = "EntMix,Main"
end

trajfile = open(args["trajfile"], "r")
if !isnothing(args["outfile"])
    outfile = open(args["outfile"], "w")
else
    outfile = open(args["trajfile"]*".stripped", "w")
end

head = read_lammpstrj_head(trajfile)

for i in 1:head["n_atoms"]+9
    readline(trajfile)
end
approx_b = position(trajfile)
secondhead = read_lammpstrj_head(trajfile)
dts = secondhead["timestep"] - head["timestep"]
maxstep = seekend_lammpstrj(trajfile)
ffstep = head["timestep"]
if args["skip"] != nothing
    @assert args["skip"] > dts "The timestep difference between the first two steps is $dts, but you specified a skipping of $(args["skip"]) which is less"
    @assert args["skip"] % dts == 0 "The timestep difference between the first two steps is $dts, but you specified a skipping of $(args["skip"]) which is not a multiple of the timestep difference"
else
    args["skip"] = dts
end

if args["startstep"] != nothing
    start = args["startstep"]
else
    seekstart(trajfile)
    start = head["timestep"]
end
if args["laststep"] != nothing
    laststep = args["laststep"]
else
    laststep,  = maxstep
end
if laststep < 0
    laststep = maxstep + laststep
end
if start < 0
    start = maxstep + start
end

find_lammpstrj_timestep(trajfile, start;delts=dts, approx_byte=approx_b)
head = read_lammpstrj_head(trajfile)
approx_b = position(trajfile)÷((head["timestep"]-ffstep)÷dts)

currentstep = head["timestep"]
@debug "Starting at step $currentstep"
writtensteps = 0
while currentstep <= laststep && !eof(trajfile)
    global approx_b
    find_lammpstrj_timestep(trajfile, currentstep;delts=dts, approx_byte=approx_b, savety=1)
    approx_b = write_step(trajfile, outfile)
    @debug "Wrote step $currentstep"
    global currentstep += args["skip"]
    global writtensteps += 1
    if args["nsteps"] != nothing && writtensteps >= args["nsteps"]
        @debug "Reached the number of steps to write, stopping"
        break
    end
end
close(trajfile)
close(outfile)
