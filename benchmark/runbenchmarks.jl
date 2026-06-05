using BenchmarkTools

include("benchmarks.jl")

function usage()
    return """
    Usage:
      julia --project=benchmark benchmark/runbenchmarks.jl [options]

    Options:
      --quick             Run each benchmark for 1 second.
      --seconds <value>   Run each benchmark for the given number of seconds.
      --output <path>     Save results as BenchmarkTools JSON.
      --list              Print benchmark names without running them.
      -h, --help          Show this message.
    """
end

function benchmark_names(group::BenchmarkGroup, prefix=String[])
    names = String[]
    for key in sort(collect(keys(group)); by=string)
        value = group[key]
        path = [prefix; string(key)]
        if value isa BenchmarkGroup
            append!(names, benchmark_names(value, path))
        else
            push!(names, join(path, " / "))
        end
    end
    return names
end

function parse_args(args)
    seconds = 5.0
    output = nothing
    list_only = false

    i = 1
    while i <= length(args)
        arg = args[i]

        if arg == "--quick"
            seconds = 1.0
        elseif arg == "--seconds"
            i == length(args) && error("Missing value after --seconds.")
            i += 1
            seconds = parse(Float64, args[i])
        elseif startswith(arg, "--seconds=")
            seconds = parse(Float64, split(arg, "=", limit=2)[2])
        elseif arg == "--output"
            i == length(args) && error("Missing value after --output.")
            i += 1
            output = args[i]
        elseif startswith(arg, "--output=")
            output = split(arg, "=", limit=2)[2]
        elseif arg == "--list"
            list_only = true
        elseif arg in ("-h", "--help")
            print(usage())
            exit(0)
        else
            error("Unknown argument: $arg")
        end

        i += 1
    end

    seconds > 0 || error("--seconds must be positive.")
    return (; seconds, output, list_only)
end

function main(args=ARGS)
    options = parse_args(args)

    if options.list_only
        println(join(benchmark_names(SUITE), "\n"))
        return nothing
    end

    results = run(SUITE; seconds=options.seconds, verbose=true)
    display(results)

    if options.output !== nothing
        BenchmarkTools.save(options.output, results)
        println("Saved benchmark results to $(options.output)")
    end

    return results
end

main()
