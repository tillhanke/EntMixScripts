# Tests for the scripts in src/, focused on calc_entropy.jl.
#
# calc_entropy.jl parses ARGS and runs main() at load time, so it is tested
# end-to-end: each test writes a small trajectory, runs the script in a fresh
# julia process and inspects its exit code, stdout/stderr and output files.
#
# Run with:  julia --project=. test/runtests.jl

using Test
using EntMix
import Chemfiles

const ROOT = dirname(@__DIR__)
const SCRIPT = joinpath(ROOT, "src", "calc_entropy.jl")
const LATTICE = "Lattice=\"10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0\" Properties=species:S:1:pos:R:3"

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

"""
Two N2 and two O2 molecules in a 10x10x10 box, shifted by `shift` along x.
Molecules are > 3 Å apart, so bond guessing finds exactly 2 molecule types.
"""
function n2o2_frame(shift=0.0)
    return [
        ("N", 1.0 + shift, 1.0, 1.0), ("N", 2.1 + shift, 1.0, 1.0),
        ("N", 6.0 + shift, 1.0, 1.0), ("N", 7.1 + shift, 1.0, 1.0),
        ("O", 1.0 + shift, 6.0, 1.0), ("O", 2.2 + shift, 6.0, 1.0),
        ("O", 6.0 + shift, 6.0, 6.0), ("O", 7.2 + shift, 6.0, 6.0),
    ]
end

"""Write `frames` (vectors of (name, x, y, z)) as xyz; with `lattice` as extended xyz."""
function write_xyz(path, frames; lattice=false)
    open(path, "w") do io
        for (i, frame) in enumerate(frames)
            println(io, length(frame))
            println(io, lattice ? LATTICE : "frame $(i - 1)")
            for (name, x, y, z) in frame
                println(io, "$name $x $y $z")
            end
        end
    end
    return path
end

"""Run calc_entropy.jl with `args`; returns (exitcode, stdout, stderr)."""
function run_script(args...; threads=1)
    out, err = IOBuffer(), IOBuffer()
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$ROOT --threads=$threads $SCRIPT $(collect(String, args))`
    proc = run(pipeline(ignorestatus(cmd); stdout=out, stderr=err))
    return proc.exitcode, String(take!(out)), String(take!(err))
end

const RESULT_LINE = r"^\s*(\d+(?:\.\d*)?),\s*(\S+)\s*$"

"""Parse "step, entropy" lines (ignoring @time output etc.) into (steps, entropies)."""
function parse_results(text)
    steps, ents = Float64[], Float64[]
    for line in split(text, '\n')
        m = match(RESULT_LINE, line)
        m === nothing && continue
        push!(steps, parse(Float64, m[1]))
        push!(ents, parse(Float64, m[2]))
    end
    return steps, ents
end

"""Run the script and return parsed results, failing the test on a nonzero exit."""
function entropies(args...; kwargs...)
    code, out, err = run_script(args...; kwargs...)
    code == 0 || @error "calc_entropy.jl failed" args err
    @test code == 0
    return parse_results(out)
end

# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

mktempdir() do dir
    plain3 = write_xyz(joinpath(dir, "plain3.xyz"), [n2o2_frame(0.5i) for i in 0:2])
    plain5 = write_xyz(joinpath(dir, "plain5.xyz"), [n2o2_frame(0.5i) for i in 0:4])
    withcell = write_xyz(joinpath(dir, "cell.xyz"), [n2o2_frame()]; lattice=true)

    @testset "calc_entropy.jl" begin

        @testset "basic run with --boxlengths" begin
            code, out, err = run_script("-b", "10,10,10", plain3)
            @test code == 0
            @test occursin("Detected 2 molecule types", err)
            steps, ents = parse_results(out)
            @test steps == [0.0, 1.0, 2.0]
            @test all(isfinite, ents)
            @test all(>(0), ents)
        end

        @testset "matches direct EntMix.entropy call" begin
            sigma = 1.3
            _, ents = entropies("-b", "10,10,10", "--sigma", string(sigma), plain3)

            traj = Chemfiles.Trajectory(plain3)
            Chemfiles.set_cell!(traj, Chemfiles.UnitCell([10.0, 10.0, 10.0]))
            first = read(traj)
            Chemfiles.type(first[0]) == "" && Molecule.type_from_name!(first)
            mols = Molecule.get_molecules(first)
            moltypes = Molecule.mol_types(first, mols)
            natoms = [length(mols[moltypes[k][1]]) for k in keys(moltypes)]
            atom_coll = [reduce(vcat, [mols[m] for m in moltypes[k]]) for k in keys(moltypes)]
            expected = [
                EntMix.entropy(Chemfiles.read_step(traj, s), atom_coll, sigma;
                               baselength="homo", dfunc=EntMix.gaus, natoms=natoms)
                for s in 0:2
            ]
            close(traj)
            @test ents ≈ expected rtol = 1e-10
        end

        @testset "step selection" begin
            @test entropies("-b", "10,10,10", plain5)[1] == [0, 1, 2, 3, 4]
            @test entropies("-b", "10,10,10", "--startstep", "2", plain5)[1] == [2, 3, 4]
            @test entropies("-b", "10,10,10", "--endstep", "1", plain5)[1] == [0, 1]
            @test entropies("-b", "10,10,10", "--stepinterval", "2", plain5)[1] == [0, 2, 4]
            @test entropies("-b", "10,10,10", "--startstep", "1", "--endstep", "3",
                            "--stepinterval", "2", plain5)[1] == [1, 3]
            # endstep past the end of the file is clamped to the last frame
            @test entropies("-b", "10,10,10", "--endstep", "100", plain3)[1] == [0, 1, 2]
        end

        @testset "multithreaded output is sorted and identical" begin
            s1, e1 = entropies("-b", "10,10,10", plain5; threads=1)
            s4, e4 = entropies("-b", "10,10,10", plain5; threads=4)
            @test s4 == s1 == [0, 1, 2, 3, 4]
            @test e4 ≈ e1 rtol = 1e-12
        end

        @testset "--outfile and --append" begin
            outfile = joinpath(dir, "out.csv")
            code, out, _ = run_script("-b", "10,10,10", "-o", outfile, plain3)
            @test code == 0
            @test isempty(parse_results(out)[1])  # results go to the file, not stdout
            steps, ents = parse_results(read(outfile, String))
            @test steps == [0, 1, 2]
            @test countlines(outfile) == 3

            # without --append the file is overwritten
            run_script("-b", "10,10,10", "-o", outfile, "--endstep", "0", plain3)
            @test parse_results(read(outfile, String))[1] == [0]

            # with --append new results are added after existing ones
            run_script("-b", "10,10,10", "-o", outfile, "--append", plain3)
            steps2, ents2 = parse_results(read(outfile, String))
            @test steps2 == [0, 0, 1, 2]
            @test ents2[2:end] ≈ ents rtol = 1e-12
        end

        @testset "unit cell handling" begin
            @testset "missing cell is an error" begin
                code, out, err = run_script(plain3)
                @test code == 2
                @test occursin("No Unit Cell", err)
                @test isempty(parse_results(out)[1])
            end

            @testset "cell read from extended xyz" begin
                from_file = entropies(withcell)[2]
                from_arg = entropies("-b", "10,10,10", write_xyz(joinpath(dir, "nocell.xyz"), [n2o2_frame()]))[2]
                @test from_file ≈ from_arg rtol = 1e-10
            end

            @testset "--boxlengths ignored when file has a cell" begin
                @test entropies("-b", "20,20,20", withcell)[2] ≈ entropies(withcell)[2] rtol = 1e-10
            end

            @testset "--boxlengths changes the result" begin
                e10 = entropies("-b", "10,10,10", plain3, "--endstep", "0")[2]
                e15 = entropies("-b", "15,15,15", plain3, "--endstep", "0")[2]
                @test !(e10 ≈ e15)
            end
        end

        @testset "--smearfunc" begin
            # sigma = 3 so that the compact-support kernels (const, linear) of
            # different molecules overlap; otherwise their entropy is exactly 0
            base = ("-b", "10,10,10", "--endstep", "0", "--sigma", "3.0")
            results = Dict(
                f => entropies(base..., "--smearfunc", f, plain3)[2][1]
                for f in ("gaus", "slater", "const", "linear")
            )
            @test all(isfinite, values(results))
            @test all(>(1e-4), values(results))
            # each kernel gives a distinct value
            @test length(unique(round.(collect(values(results)); sigdigits=6))) == 4
            # default kernel is gaus
            @test entropies(base..., plain3)[2][1] ≈ results["gaus"]

            code, out, err = run_script("-b", "10,10,10", "--smearfunc", "nonsense", plain3)
            @test code == 2
            @test occursin("nonsense", err)
            @test isempty(parse_results(out)[1])
        end

        @testset "--sigma and --sigmatype" begin
            base = ("-b", "10,10,10", "--endstep", "0")
            homo = entropies(base..., plain3)[2][1]
            @test entropies(base..., "--sigmatype", "homo", plain3)[2][1] ≈ homo
            vdw = entropies(base..., "--sigmatype", "VDW", plain3)[2][1]
            cov = entropies(base..., "--sigmatype", "Covalent", plain3)[2][1]
            @test isfinite(vdw) && isfinite(cov)
            @test !(vdw ≈ homo) && !(cov ≈ homo) && !(vdw ≈ cov)

            wide = entropies(base..., "--sigma", "2.0", plain3)[2][1]
            @test !(wide ≈ homo)
        end

        @testset "mixed configuration has higher entropy than segregated" begin
            # same molecule positions, only the species assignment differs
            sites = [(2.0, 2.0, 2.0), (2.0, 5.0, 2.0), (7.0, 2.0, 7.0), (7.0, 5.0, 7.0)]
            mol(el, (x, y, z)) = [(el, x, y, z), (el, x + (el == "N" ? 1.1 : 1.2), y, z)]
            mixed = reduce(vcat, mol.(["N", "O", "O", "N"], sites))
            segregated = reduce(vcat, mol.(["N", "N", "O", "O"], sites))
            s_mix = entropies("--sigma", "1.5", write_xyz(joinpath(dir, "mix.xyz"), [mixed]; lattice=true))[2][1]
            s_seg = entropies("--sigma", "1.5", write_xyz(joinpath(dir, "seg.xyz"), [segregated]; lattice=true))[2][1]
            @test s_mix > s_seg > 0
        end

        @testset "argument errors" begin
            # missing positional trajfile
            code, _, err = run_script()
            @test code != 0
            @test occursin("trajfile", err)
            # nonexistent trajectory
            code, out, _ = run_script("-b", "10,10,10", joinpath(dir, "does_not_exist.xyz"))
            @test code != 0
            @test isempty(parse_results(out)[1])
        end
    end
end
