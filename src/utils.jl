function save1D(storage::Array, savename::String, path::String; headername="data")
	n = length(storage)
	t=0
	for field in storage
		CSV.write(datadir(path*savename*"t=$t"), DataFrame(headername = field), writeheader=false)
	t += 1
    end
end
