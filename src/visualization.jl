
function plotAnimation(storage::Array, gifname::String)
    t = 0
	n = length(storage)
	p = ProgressMeter.Progress(n, 0.5, "Creating a gif...")
	anim = @animate for field in storage
        plot(field, ylim = (0, 1), label = "time=$t")
		ProgressMeter.next!(p)
        t += 1
    end

    gif(anim, gifname, fps = 15)
end