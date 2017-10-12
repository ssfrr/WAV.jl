# FileIO integration support
loadstreaming(s::Stream{format"WAV"}) = WAVSource(s.io)
load(s::Stream{format"WAV"}) = read(WAVSource(s.io))
savestreaming(s::Stream{format"WAV"}; kwargs...) = WAVSink(s.io, kwargs...)
save(s::Stream{format"WAV"}, data; kwargs...) = write(WAVSink(io, kwargs...), data)

# TODO: think about who's responsible for closing this file
loadstreaming(f::File{format"WAV"}; kwargs...) = WAVSource(open(f.filename))
function load(f::File{format"WAV"}; kwargs...)
    open(f.filename) do io
        read(WAVSource(io))
    end
end

# TODO: think about who's responsible for closing this file
savestreaming(f::File{format"WAV"}; kwargs...) = WAVSink(open(f.filename, "w"))
function save(f::File{format"WAV"}, data; kwargs...)
    open(f.filename, "w") do io
        write(WAVSink(io, kwargs...), data)
    end
end
