
# This will take a text file of space-separated RGB values and create the
# proper `register_gradient_colors(:magma, sample_evenly([...]))` calls

function build_grad_string(fn::AbstractString, libname, separator = ' ')
    name = splitext(basename(fn))[1]
    if length(name) > 4 && name[end-3:end] == "-rgb"
        name = name[1:end-4]
    end
    @show name
    colors = [begin
        "    RGB($(join(split(strip(l), separator), ", "))),\n"
    end for l in eachline(open(fn))]
    "\nregister_gradient_colors(:$(name), sample_evenly([\n$(join(colors))], $(min(64, length(colors)))), :$(libname))\n"
end

function build_grad_file(dir::AbstractString, outname::AbstractString, separator = ' ')
    outfn = joinpath(dirname(@__FILE__()), outname*".jl")
    @show outfn
    outfile = open(outfn, "w")
    write(outfile, "register_color_library(:$(outname), ColorLibrary(Dict(:default => :sequential)))\n")
    for (root, dirs, files) in walkdir(dir)
        @show root, dirs, files
        for fn in files
            #if length(fn) > 8 && fn[end-7:end] == "-rgb.txt" #specific to cmocean
            if fn[end-3:end] == ".csv"
                write(outfile, build_grad_string(joinpath(root, fn), outname, separator))
            end
        end
    end
    close(outfile)
end
