# We use [dcraw](http://www.cybercom.net/~dcoffin/dcraw/) to extract from RAW images
# the blue channel pixels and convert to 16bit pgm without exposure manipulation. 

# The website notes: "Unless otherwise noted in the source code, these programs 
# are free for all uses, although I would like to receive credit for them."
# The program dcraw is written by Dave Coffin.


# On mac, compile dcraw.c with
# 	`llvm-gcc -o dcraw dcraw.c -lm -DNO_JPEG -DNO_LCMS -DNO_JASPER`
# source [http://vkphotoblog.blogspot.be/2014/05/dcraw-921-for-os-x-mavericks-users.html]

# ### dcraw options
# * -W: avoids stretching ("No matter how dark an image is, dcraw's auto-exposure stretches it so 
#       that one percent of its pixels appear white. The "-W" option avoids this behavior.")
# * -4: is for linear 16-bit, same as -6 -W -g 1 1 (with -g for gamma correction)
# * -j: don't stretch or rotate raw pixels
# * -t [0-7]: flip image (0=none)
# * -D: document mode without scaling (totally raw), while -d scales to 16bit (eg from 14bit). 
#       Either use -D and then reinterpret in julia or use -d. We now use -d otherwise we need to 
#       extract image property for bit depth (12bit, 14bit, ...).
# * -r 0 0 0 1: select only the blue channel (see Brusa and Bunker, 2014). This option selects from a bayer RGGB layout.
# * -v: for verbose output


const OVEREXP = 0.005

const RAW_EXT = String[".3fr", ".ari", ".arw", ".bay", ".crw", ".cr2",
".cap", ".dcs", ".dcr", ".dng",
".drf", ".eip", ".erf", ".fff", ".iiq", ".k25", ".kdc", ".mdc", ".mef", ".mos", ".mrw",
".nef", ".nrw", ".obm", ".orf", ".pef", ".ptx", ".pxn", ".r3d", ".raf", ".raw", ".rwl",
".rw2", ".rwz", ".sr2", ".srf", ".srw", ".tif", ".x3f"]

const DCRAW_DIR = joinpath(dirname(pathof(LeafAreaIndex)), "dcraw")

@static if Sys.iswindows() 
    const DCRAW_EXE = joinpath(DCRAW_DIR, "dcraw-9.26-ms-64-bit.exe")
end

@static if Sys.isunix()
    const DCRAW_EXE = joinpath(DCRAW_DIR, "dcraw")
    is_user_executable(file) = isodd(uperm(file))
    is_user_executable(DCRAW_EXE) || chmod(DCRAW_EXE, 0o755)
end

@assert isfile(DCRAW_EXE)

"""Read in the raw blue channel from a raw image. The image is copied to
a temporary directory and run through dcraw with options `-d -4 -j -t 0 -r 0 0 0 1`."""
rawblueread(filepath::String; kwargs...) = readraw(filepath; kwargs...)

rawcolourread(filepath::String; kwargs...) = readraw(filepath; colour=true, kwargs...)

function readraw(filepath::String; overwrite=false, rmcopy=true, rmpxm=true,
    destdir=joinpath(tempdir(), "pxm"), colour=false, kwargs...)

    @assert isfile(filepath) 
    ext = lowercase(last(splitext(filepath)))
    @assert ext ∈ RAW_EXT
    isdir(destdir) || mkdir(destdir)

    if colour
        raw2pxm = f -> run(`$DCRAW_EXE -4 -j -t 0 $f`)
        pxm_ext = ".ppm"
    else
        #raw2pxm = f -> run(`$DCRAW_EXE -d -4 -j -t 0 -r 0 0 0 1 $f`)
        #pxm_ext = ".pgm"
		# add De-Bayering
		raw2pxm = f -> run(`$DCRAW_EXE -4 -j -t 0 -q 3 -o 0 $f`)
		pxm_ext = ".ppm"
    end

    copyfile = joinpath(destdir, last(splitdir(filepath)))
    pxm = first(splitext(copyfile)) * pxm_ext

    if !overwrite && isfile(pxm)
        return FileIO.load(pxm)
    end

    cp(filepath, copyfile; force=true)
    raw2pxm(copyfile)
    img = FileIO.load(pxm)
    rmcopy && rm(copyfile)
    rmpxm && rm(pxm)
    img
end

isoverexposed(img::AbstractArray) = sum(img .== 1) > OVEREXP * length(img)
isoverexposed(polim::PolarImage) = isoverexposed(pixels(polim))