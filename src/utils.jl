"""
    save1D(storage::Array, savename::String, path::String; headername="data")

writes a one dimensional field to disk by providing the field in `storage`, a savename as `String` and a path.
"""
function save1D(storage::Array, savename::String, path::String; headername="data")
    n = length(storage)
    t=0
    for field in storage
        CSV.write(datadir(path*savename*"t=$t"), DataFrame(headername = field), writeheader=false)
        t += 1
    end
end
