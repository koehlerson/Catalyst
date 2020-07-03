
"""
	plotAnimation(storage::Array, gifname::String, lim=(0,1))

Plots a collection of one dimensional arrays to a gif. Applicable to plot the
Reaction Operator or the one dimensional concentration field at each time instance 
over the domain. Both mentioned things are returned by Catalyst.solve(...) 
"""
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

"""
	plotOverTime(storage::Array, ylim=(0,1), save=false, figname=missing)

takes the last entry of the `storage` array in each time step and plots it over time.
Useful to plot the concentration field returned by Catalyst.solve() over time.
"""
function plotOverTime(storage::Array, ylim=(0,1), save=false, figname=missing)
	n = length(storage[1])
	plot(getindex.(storage, n), ylims=ylim)
end
