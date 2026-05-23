using Potential
using Test
using Logging

@testset "Potential.jl" begin
    @testset "read_particle_coords" begin
        mktempdir() do dir
            expected = [0.0 1.0 2.0; 3.0 4.0 5.0]

            whitespace_path = joinpath(dir, "coords.txt")
            write(whitespace_path, "# x y z\n0 1 2\n3\t4\t5\n")
            @test read_particle_coords(whitespace_path) == expected
            @test read_particle_coords(whitespace_path; factor=2) == 2 .* expected

            csv_path = joinpath(dir, "coords.csv")
            write(csv_path, "x,y,z\n0,1,2\n3,4,5\n")
            @test read_particle_coords(csv_path) == expected

            tsv_path = joinpath(dir, "coords.tsv")
            write(tsv_path, "0\t1\t2\n3\t4\t5\n")
            @test read_particle_coords(tsv_path) == expected
        end
    end

    @testset "basic numerical helpers" begin
        m, c = linear_fit([1.0, 2.0, 3.0], [3.0, 5.0, 7.0])
        @test m ≈ 2.0
        @test c ≈ 1.0

        coords = [0.0 0.0 0.0; 1.0 0.0 0.0]
        sf = structure_factor_debye_histogram(coords; qmin=0.1, qmax=0.2, nq=3, show_progress=false)
        @test length(sf.q) == 3
        @test length(sf.Sq) == 3
        @test all(isfinite, sf.Sq)
    end

    @testset "threaded RDF" begin
        if Threads.nthreads() > 1
            coords = [0.0 0.0 0.0;
                      1.0 0.0 0.0;
                      0.0 1.0 0.0;
                      0.0 0.0 1.0;
                      1.0 1.0 1.0]
            r = range(0.0, stop=2.0, length=12)
            @test rdf_analytic_norm(coords, r, 0.2; threaded=true) ≈
                  rdf_analytic_norm(coords, r, 0.2; threaded=false)
            @test rdf_analytic_norm(coords[1:3, :], coords[3:5, :], r, 0.2; threaded=true) ≈
                  rdf_analytic_norm(coords[1:3, :], coords[3:5, :], r, 0.2; threaded=false)
        else
            @test true
        end
    end

    @testset "Julia version support warning" begin
        @test Potential._is_tested_julia_version(v"1.12.0")
        @test Potential._is_tested_julia_version(v"1.12.99")
        @test !Potential._is_tested_julia_version(v"1.11.0")
        @test !Potential._is_tested_julia_version(v"1.13.0")

        msg = Potential._untested_julia_warning_message(v"1.13.0")
        @test occursin("Potential.jl has not been tested on Julia 1.13", msg)
        @test occursin("tested on Julia 1.12", msg)

        @test_logs (:warn, r"Potential\.jl has not been tested on Julia 1\.13") Potential._warn_if_untested_julia_version(v"1.13.0")
        @test_logs min_level=Logging.Warn Potential._warn_if_untested_julia_version(v"1.12.0")
    end
end
