using Documenter, LeafAreaIndex

makedocs(
    modules = [LeafAreaIndex],
    clean    = true,
    format   = :html,
    sitename = "Leaf Area Index Documentation",
    authors  = "Ken Bastiaensen",
    pages    = [
        "LeafAreaIndex.jl" => "index.md",
        "Quick Intro"      => "quick_intro.md",
        "Gap Fraction"     => "gapfraction.md",
        "LAI Inversion"    => "LAI.md",
        "Clumping"         => "clumping.md",
        "Calibration"      => "calibration.md",
    ],
    linkcheck = false,#!("skiplinks" in ARGS),
    # Use clean URLs, unless built as a "local" build
    html_prettyurls = false#!("local" in ARGS)
)
deploydocs(
    #root   = Pkg.dir("LeafAreaIndex","docs"),
    repo   = "github.com/ETC-UA/LeafAreaIndex.jl.git",
    target = "build",
    julia  = "0.6", 
    deps   = nothing, 
    make   = nothing
)
