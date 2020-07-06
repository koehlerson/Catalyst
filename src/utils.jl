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

@doc raw"""
    save_all(c::Array, R::Array, savename::String, path::String)

stores the one dimensional scalar fields `c` and `R` as `"c_fields_$savename"` in the directory
`"$path/c/"`, `"$path/R/"` and `"$path/cPsim/"`, respectively 
"""
function save_all(c::Array, R::Array, savename::String, path::String)
    save1D(c, "c_field_$savename", "$path/c/")
    save1D(R, "R_field_$savename", "$path/R/")
    number_of_nodes = length(c[1])
    cPsim = getindex.(c, number_of_nodes)
    CSV.write("$path/cPsim/cPsim_$savename", DataFrame(headername=cPsim), writeheader=false)
end
