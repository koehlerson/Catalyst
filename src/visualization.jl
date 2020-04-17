
function plotAnimation(storage::Array, gifname::String, lim=(0,1))
    t = 0
	n = length(storage)
	p = ProgressMeter.Progress(n, 0.5, "Creating a gif...")
	anim = @animate for field in storage
        plot(field, ylim = lim, label = "time=$t")
		ProgressMeter.next!(p)
        t += 1
    end

    gif(anim, gifname, fps = 24)
end

function plotOverTime(storage::Array, ylim=(0,1), save=false, figname=missing)
		n = length(storage[1])
		plot(getindex.(storage, n), ylims=ylim)
end
