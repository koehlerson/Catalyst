
function plotAnimation(storage::Array, gifname::String)
    anim = @animate for c in store
        plot(c)
    end

    gif(anim, gifname, fps = 30)
end

export plotAnimation
