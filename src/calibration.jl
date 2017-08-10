
# some auxilary functions
dist(p::Point, q::Point) = norm(p-q)
angle(p::Point, q::Point) = atan2(p.y - q.y, q.x - p.x)

# Center
# ------

# Calibrates the lens center for a camera+lens setup. Each combination of a lens and camera
# is assumed to be unique and every lens center should be calibrated.

# Luckily this is a straightforward procedure, taken from the
#[CAN-EYE user guide](http://www6.paca.inra.fr/can-eye/Documentation-Publications/Documentation).
# Dril 3-4 small holes in your lens cap (or prick an aluminium foil well attached to a lens
#ring) and take a picture for each rotation. Because coordinates of the dots in the pictures
# should lie on circle, the lens center can be calculated.

# The dataframe requires 3 columns with resp. the x-coordinates, corresponding y-coordinates
# and a tag (eg a number) for each dot/circle.
function calibrate_center(df::DataFrame, height::Int, width::Int)

    # First we check the dataframe input and create Nρ and Nϕ, resp. number of circles and
    # points per circle.
    @assert width >= height
    for name = [:x, :y, :circle]
        name in names(df) || error("Column $name not found in dataframe input.")
    end
    eltype(df[:x]) <: Real || throw(DomainError("Column :x should be numeric."))
    eltype(df[:y]) <: Real || throw(DomainError("Column :y should be numeric."))

    pool!(df, :circle)
    Nρ = length(levels(df[:circle]))
    @assert 3 <= Nρ <= 5

    all(by(df, :circle, length)[2]) > 0 || error("Different number of pictures per circle.")
    Nϕ = Int(nrow(df) / Nρ)
    @assert Nϕ * Nρ == nrow(df)

    @assert all(df[:x] .< width)
    @assert all(df[:y] .< height)

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
    pic_center = Point(floor(Int, width/2),floor(Int, height/2))

    # We optimize using the Nelder-Mead downhill simplex algorithm,
    # the default for functions without a derivate explicetly defined.
    res = optimize(sse, initsse(pic_center, points), iterations=10_000)
    return lensx, lensy = round(Int, minimizer(res[1:2]))
end


# Projection Function
# -------------------

const N_LATERAL = 5 # number of points on side taken for Δ estimation

# TODO use lensx, lensy from calibrated lens center?
function calibrate_projfun(df::DataFrame, height::Int, width::Int)

    # First we check the dataframe inputs
    @assert width >= height 
    for name = [:cm, :px, :H, :pos]
        @assert name in names(df)
    end
    eltype(df[:cm]) <: Real || throw(DomainError("Column :cm should be numeric."))
    eltype(df[:px]) <: Real || throw(DomainError("Column :px should be numeric."))

    @assert all(df[:px] .< width)

    #pool!(df, [:H, :pos])
    @assert length(levels(df[:H])) == 2
    @assert length(levels(df[:pos])) == 3
    for name = ["right", "left", "perpendicular"]
        @assert name in levels(df[:pos])
    end


    # We use 2 pictures at different distances H1, H2 from the perpendicular ruler to
    # estimate the distance to the perpendicular ruler more precisely than with a physical
    # ruler during setup.

    # We first estimate more correctly the distance difference Δ between im1 and im2 by
    # noticing that the shift on the lateral rulers between the two pictures is equal to Δ.
    @show H1, H2 = extrema(df[:H])

    c1 = df[(df[:pos] .== "perpendicular") & (df[:cm] .== 0) & (df[:H] .== H1), :px][1]
    c2 = df[(df[:pos] .== "perpendicular") & (df[:cm] .== 0) & (df[:H] .== H2), :px][1]

    # Utility function to interpolate a pixel from image 1 to the ruler of image 2
    function interpol(x, df2, center2)
        # suffix 2 is for second image with H2
        s2 = sort(df2, cols=:px)
        x2 = s2[:px] - center2
        cm2 = s2[:cm]

        if (x < minimum(x2)) | (x > maximum(x2))
            return NaN
        end

        ind = searchsortedfirst(x2, x)
        out = cm2[ind] - (x2[ind] - x) / (x2[ind] - x2[ind-1]) * (cm2[ind] - cm2[ind-1])
    end

    # fitting function to optimize Δ
    function fitΔ(Δ, df1, df2, center1, center2; errfun::Function=x->x^2)
        diffs = Float64[]
        for row in eachrow(df1)
            x = row[:px] - center1
            cm2 = interpol(x, df2, center2)
            !isnan(cm2) && push!(diffs, cm2 - row[:cm] - Δ)
        end
        return norm(diffs)
    end

    # starting point for optimization
    Δinit = H2 - H1

    left1 = head(sort!(df[(df[:pos] .== "left") & (df[:H] .== H1), :], cols=:cm, rev=true),
                N_LATERAL)
    left2 = head(sort!(df[(df[:pos] .== "left") & (df[:H] .== H2), :], cols=:cm, rev=true),
                N_LATERAL)
    resleft = optimize(Δ -> fitΔ(Δ, left1, left2, c1, c2), Δinit/2, Δinit * 3/2)

    right1 = head(sort!(df[(df[:pos] .== "right") & (df[:H] .== H1), :], cols=:cm, rev=true),
                N_LATERAL)
    right2 = head(sort!(df[(df[:pos] .== "right") & (df[:H] .== H2), :], cols=:cm, rev=true),
                N_LATERAL)
    resright = optimize(Δ -> fitΔ(Δ, right1, right2, c1, c2), Δinit/2, Δinit * 3/2)

    Δ = (minimizer(resleft) + minimizer(resright)) / 2
    @show Δ

    # With a more precise estimation of Δ, we can estimate H1 more precisely.
    # We optimize H1 for the value that provides the best fit between the two projection
    # functions of the two images. We fit a 2nd degree polynomial without intercept to the
    # projection functions ρ(θ) = aθ + bθ^2.
    # For higher order polynomial, simply adjust the `abfit` function accordingly.

    # perpendicular ruler half width
    w = maximum(df[df[:pos] .== "perpendicular", :cm ])

    function projfunpoints(df, center, H)
        θset = Float64[]
        ρset = Int[]
        for row in eachrow(df)
            cm = row[:cm]
            ρ = abs(row[:px] - center)
            θ = row[:pos] == "perpendicular" ? atan(abs(cm) / H) : atan(w / (H - cm))
            push!(θset, θ)
            push!(ρset, ρ)
        end
        θset, ρset
    end

    function abfit(θ, ρ)
    # fit ρ(θ) = aθ + bθ²
        A = [θᵢ^p for θᵢ in θ, p = 1:2]
        a, b = A \ ρ
    end

    # plot to spot outliers
    #function plotfit(H,c)
    #    lens = df[df[:H] .== H, :]
    #    θ, ρ = projfunpoints(lens, c, H)
    #    Winston.plot(θ, ρ, "o")
    #    a, b = abfit(θ, ρ)
    #    Winston.oplot(θ -> a*θ + b*θ.^2, 0, pi/2, "r")
    #    Winston.title("Observations and fit for $H")
    #end
    #display(plotfit(H1, c1)); display(plotfit(H2, c2))

    function projfundiff(H, Δ, lens1, lens2, center1, center2)
        a1, b1 = abfit(projfunpoints(lens1, center1, H)...)
        a2, b2 = abfit(projfunpoints(lens2, center2, H+Δ)...)
        (a1 - a2).^2 + (b1-b2).^2
    end

    lens1 = df[df[:H] .== H1, :]
    lens2 = df[df[:H] .== H2, :]
    Hres = optimize(H->projfundiff(H, Δ, lens1, lens2, c1, c2), H1 / 2, H1 * 2)
    H = minimizer(Hres)
    @show H

    a1, b1 = abfit(projfunpoints(lens1, c1, H)...)
    a2, b2 = abfit(projfunpoints(lens2, c2, H+Δ)...)

    a, b = (a1 + a2) / 2, (b1 + b2) / 2
    return @show a, b
end
