using Base.Graphics: Point, norm
using Optim
using DataFrames
# some auxilary functions
dist(p::Point, q::Point) = norm(p-q)
angle(p::Point, q::Point) = atan2(p.y - q.y, q.x - p.x)

function calibrate_center(df::DataFrame, height::Int, width::Int)

    # First we check the dataframe input and create Nρ and Nϕ, resp. number of circles and
    # points per circle.
    width >= height || error("Width is smaller than height.")
    for name = [:x, :y, :circle]
        name in names(df) || error("Column $name not found in dataframe input.")
    end
    eltype(df[:x]) <: Real || throw(DomainError("Column :x should be numeric."))
    eltype(df[:y]) <: Real || throw(DomainError("Column :y should be numeric."))

    pool!(df, :circle)
    Nρ = length(levels(df[:circle]))
    3 <= Nρ <= 5 || error("3 to 5 dots per picture expected.")

    all(by(df, :circle, length)[2]) > 0 || error("Different number of pictures per circle.")
    Nϕ = int(nrow(df) / Nρ)
    Nϕ * Nρ == nrow(df) || error("Inconsistent dataframe row count.")

    all(df[:x] .< width)  || error("Some x coordinates larger than width.")
    all(df[:y] .< height) || error("Some y coordinates larger than height.")

    # Create Matrix of Points with a column per circle.
    points = Point[]
    by(df, :circle) do subdf
        for row in eachrow(subdf)
            push!(points, Point(row[:x], row[:y]))
        end
    end
    points = reshape(points, (Nϕ, Nρ))

    # Create fit function as sum of squared errors.
    function sse(inp::Vector)
        # all parameters in single inp vector
        cx, cy = inp[1:2]
        c = Point(cx, cy)
        ρ = inp[3:(Nρ + 2)] #circle radius
        δ = inp[(3 + Nρ):(2Nρ + 2)] #angle offset between circles in single picture
        ϕ = inp[(3 + 2Nρ):end] #angle offset between pictures

        se = Float64[]
        for i = 1:Nρ
            # Normalize this squared error with squared typical radius because we sum upe
            # the errors for both radius and angles.
            push!(se, sum([(ρ[i] - dist(c,points[n, i]))^2 for n = 1:Nϕ])/(height^2/4))
        end
        for i = 1:Nϕ
            push!(se, sum([(ϕ[i] + δ[n] - angle(c,points[i, n]))^2 for n = 1:Nρ]))
        end
        sum(se)
    end

    # Estimate the starting point of the opimization from the inner circle.
    function initsse(c::Point, points::Matrix{Point})
        init = zeros(2 + 2Nρ + Nϕ)
        init[1:2] = [c.x, c.y]
        for i = 1:Nρ
            init[2+i] = dist(c, points[1,i])
            #offset
            init[2+Nρ+i] = angle(c, points[1,i]) - angle(c,points[1,1])
        end
        for i = 1:Nϕ
            init[2+2Nρ+i] = angle(c, points[i,1])
        end
        init
    end
    pic_center = Point(ifloor(width/2),ifloor(height/2))

    # We optimize using the Nelder-Mead downhill simplex algorithm,
    # the default for functions without a derivate explicetly defined.
    res = optimize(sse, initsse(pic_center, points), iterations=10_000)
    return lensx, lensy = int(res.minimum[1:2])
end
