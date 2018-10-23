var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "LeafAreaIndex.jl",
    "title": "LeafAreaIndex.jl",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#LeafAreaIndex.jl-1",
    "page": "LeafAreaIndex.jl",
    "title": "LeafAreaIndex.jl",
    "category": "section",
    "text": "The LeafAreaIndex.jl package offers several methods to work with hemispherical pictures for the determination of Leaf Area Index (LAI).The ICOS - ETC team at the University of Antwerp provides a free online service with automatic LAI estimation: https://icos.ua.ac.be/For more information on the organization see the website of ICOS Belgium.The code (and this documentation) is licensed under the permissive MIT license."
},

{
    "location": "quick_intro.html#",
    "page": "Quick Intro",
    "title": "Quick Intro",
    "category": "page",
    "text": ""
},

{
    "location": "quick_intro.html#LAI-1",
    "page": "Quick Intro",
    "title": "LAI",
    "category": "section",
    "text": "For an introduction on LAI we refer to Weiss et al., 2004[1] and Thimonier et al., 2010[2]."
},

{
    "location": "quick_intro.html#Quick-Introduction-1",
    "page": "Quick Intro",
    "title": "Quick Introduction",
    "category": "section",
    "text": "Install the package throughPkg.clone(\"https://github.com/ETC-UA/LeafAreaIndex.jl\")You construct a PolarImage from a Calibration type and an Image (or in general, a Matrix). The calibration requires the image size, the coordinates of the lens center and the (inverse) projection function.  (The projection function maps polar distance ρ [in pixels] on the image to the zenith angle θ [in radians] of the scene and is usually not linear.)using Images\nimg = imread(\"image.dng\")\nimgblue = blue(img) #take the blue channel\nusing LeafAreaIndex\nmycameralens = calibrate(height, width, ci, cj, fθρ, fρθ)\npolarimg = PolarImage(imgblue, mycameralens)The first step is automatical thresholding with the default method Ridler Calvard:thresh = threshold(polarimg)The clumping factor can be estimated with the method of Lang Xiang with 45ᵒ segments between view angles θ1 and θ2 using langxiang45(polarimg, thresh, θ1, θ2). Similarly chencihlar(polarimg, thresh, θ1, θ2) for the Chen Chilar method.There are two methods to estimate LAI assuming an ellipsoidal leave angle distribution, both also estimating the Average Leaf Inclination Angle (ALIA). The first one uses a Lookup Table (LUT) and the second one an optimization method.LAI1 = ellips_LUT(polarimg, thresh)\nLAI2 = ellips_opt(polarimg, thresh)"
},

{
    "location": "quick_intro.html#further-methods-1",
    "page": "Quick Intro",
    "title": "further methods",
    "category": "section",
    "text": "For images taken (vertically upwards) on a domain with slope of eg 10ᵒ and downward to the East, you must include this in your PolarImage.myslope = Slope(10/180*pi, pi/2)\npolarimg = PolarImage(img, mycameralens, myslope)Two more automatic thresholding methods can be used with edge_threshold and minimum_threshhold for the Edge Detection and Minimum algorithm.Further LAI estimation methods are available: zenith57(polarimg, thresh) for using a view angle of 57ᵒ where the ALIA influence is minimal; miller(polarimg, thresh for miller integration assuming a constant leaf angle; and lang(polarimg, thresh) that uses a first order regression on the miller method.To access the pixels in a particular zenith range, pixels(polarimg, pi/6, pi/3) will return a vector with pixels quickly, sorted by increasing ρ (and then by polar angles ϕ for identical ρ). A shortcut pixels(polarimg) is translated to pixels(polarimg, 0, pi/2).The segments function can further split these ring pixels in n segments (eg. for clumping calculation). It returns a vector with n elements, each a vector with segment pixels.You can also construct an iterator to access a specific zenith range. It will return the pixels on each ring in the range by increasing integer ρ² in a tuple with a vector of polar angles ϕ and a vector of corresponding pixels.for (ρ², ϕ, px) in rings(polarimg, pi/6, pi/3)\n    # do something with each ρ², ϕ, px variable\nendIn case of problems or suggestion, don\'t hesitate to submit an issue through the issue tracker or code suggestions through a pull request.[1]: Weiss et al 2004[2]: Thimonier et al 2010"
},

{
    "location": "gapfraction.html#",
    "page": "Gap Fraction",
    "title": "Gap Fraction",
    "category": "page",
    "text": ""
},

{
    "location": "gapfraction.html#Gap-Fraction-1",
    "page": "Gap Fraction",
    "title": "Gap Fraction",
    "category": "section",
    "text": ""
},

{
    "location": "gapfraction.html#Binarization-1",
    "page": "Gap Fraction",
    "title": "Binarization",
    "category": "section",
    "text": "In a first step each pixel of the image is categorized as either sky or vegetation. This is called binarization. The binarization methods in this package all require the ability to order image pixels (technically, to be able to compare two pixels and determine if one\'s values is less than the other\'s). Therefore we work on the basis of a single light intensity value for each pixel, as this is not defined unambiguously for RGB values. We recommend using the Image package to select the  blue channel using the blue method. If you wish, however, you can easily use any linear combination of the different colour channels.In general, these binarization methods determine a single global threshold between sky and vegatation. Any pixel with values above this threshold is considered sky.The standard method is from Ridler Calvard, 1978, which is similar to k-means clustering with two classes in one dimension.The second method is the Minimum method that looks for the minimum between two peaks in the histogram.Finally, the Edge Detection method by Nobis & Hunziker, 2004 aims to maximize the contrast at edges between sky and vegetation.A comparison of different methods can be found in Glatthorn & Beckschäfer, 2014.You can access the different methods with (where polimg is a PolarImage instance):thresh_RC = threshold(polimg) or thresh_RC = RidlerCalvard(polimg)\nthresh_min = minimum_threshold(polimg)\nthresh_edge = edge_threshold(polimg)The RidlerCalvard and minimum_threshold only use pixels within a 90ᵒ view angle, while edge_threshold cuts out a rectangle around the visible pixels to maintain order of neighbouring pixels."
},

{
    "location": "gapfraction.html#Gap-Fraction-2",
    "page": "Gap Fraction",
    "title": "Gap Fraction",
    "category": "section",
    "text": "Once a threshold ahs been set, you can calculate a gap fraction for a specific part of the image and specialized methods exist for its logarithm and the contact frequency. However, these low level functions are in general not required to be called by the user but are automatically called by higher level functions. See the implementation section for more details."
},

{
    "location": "LAI.html#",
    "page": "LAI Inversion",
    "title": "LAI Inversion",
    "category": "page",
    "text": ""
},

{
    "location": "LAI.html#LAI-estimation-1",
    "page": "LAI Inversion",
    "title": "LAI estimation",
    "category": "section",
    "text": "This package provides several methods to estimate LAI. All of them are based on the same inversion model."
},

{
    "location": "LAI.html#Model-1",
    "page": "LAI Inversion",
    "title": "Model",
    "category": "section",
    "text": "We estimate LAI from the gap fraction through inversion of the Poisson model:T(theta_V) = expbig(-k(theta_V) L_ebig) = expBig(-fracG(theta_V)L_ecostheta_VBig)with T(theta_V): the observed transmission or gap fraction,\ntheta_V: the view zenit angle,\nL_e: the effective Leaf Area Index with L_e = L cdotOmega and Omega the cluming factor/correction,\nk(theta_V): the extinction coefficient that depends on G(theta_V), and\nG(theta_V): the projection function, which in its turn depends on a leaf angle distribution.Further, the contact frequency K is defined as K(theta_V) = G(theta_V) cdot L_e = -log T(theta_V) cos theta_V and is used as observed value to compare with the modelled value for LAI fitting."
},

{
    "location": "LAI.html#Projection-function-1",
    "page": "LAI Inversion",
    "title": "Projection function",
    "category": "section",
    "text": "The projection function accounts for the interaction between the incoming beam of light from direction (theta_V phi_V) with the leaves assuming a certain leaf angle distribution with leave angles (theta_L phi_L). More specifically, the projection function is the mean projection of a unit foliage area in the direction (theta_V phi_V).In general the projection function G(theta_V phi_V) is definded as:G(theta_V phi_V) = frac12piint_0^2piint_0^pi2lvertcospsirvert g(theta_Lphi_L)sintheta_L mathrmdtheta_Lmathrmdphi_L with cospsi = costheta_Vcostheta_L + sintheta_Vsintheta_Lcos(phi_V - phi_L).For forest canopies the azimuth leaf angle phi_L distribution is often assumed constant which reduces G toG(theta_V) = int_0^pi2A(theta_Vtheta_L)g(theta_L)mathrmdtheta_Lwith a direction function A: A(theta_V theta_L) = begincases \n    costheta_Vcostheta_L  quadlvertcottheta_Vcottheta_Lrvert1  \n    costheta_Vcostheta_L1+tfrac2pi(tanbeta - beta) quadtextotherwise\nendcaseswith beta = cos^-1(cottheta_Vcottheta_L).The projection function G(theta_V) as well as the leaf angle distribution function g(theta_L) must obey the following normalization functions:int_0^pi2G(theta_V)sintheta_Vmathrmdtheta_V = tfrac12\nint_0^pi2g(theta_V)sintheta_Vmathrmdtheta_V = 1"
},

{
    "location": "LAI.html#Leaf-Angle-Distribution-1",
    "page": "LAI Inversion",
    "title": "Leaf Angle Distribution",
    "category": "section",
    "text": "The standard leaf angle distribution is the _ellipsoidal_ distribution with 1 parameter chi:g(theta_L chi) = frac2chi^3sintheta_LD(cos^2theta_L + chi^2sin^2theta_L)^2 with D approx chi + 1774(chi+1182)^-0733.The resulting projection function then becomes (formula(A.6) in Thimonier et al. 2010[2]):G(theta_V chi) = fraccostheta_Vsqrtchi^2 + tan^2theta_Vchi + 1702(chi+112)^-0708Thimonier et al. 2010[2]: If the vertical semi-axis is a and the horizontal semi-axis b, the ellipsoidal leaf angle distribution parameter is defined as chi = b  a. The parameter chi is directly related to the average leaf inclination angle (ALIA) bartheta_L through formula (30) in Wang et al 2007:chi approx Big(fracbartheta_L965Big)^-06061 - 3As an example we plot the ellipsoidal for different average leaf angles:g(θL, χ) = 2χ^3 * sin(θL) /( (χ + 1.774(χ + 1.182)^-0.733) * (cos(θL)^2 +  χ^2*sin(θL)^2)^2)\nALIA_to_x(ALIA) = (ALIA/9.65).^-0.6061 - 3\nusing Winston\np = plot()\nfor alia in linspace(5 *pi/180, 85 *pi/180, 10)\n    oplot(θL -> g(θL, ALIA_to_x(alia)), 0, pi/2, \"--\")\nend\ntitle(\"leaf angle distribution \\\\theta_L for different average leaf angles\")\nylim(0, 5) ; xlabel(\"leaf angle \\\\theta_L\")(Image: ellipsoidal leaf angle distribution)"
},

{
    "location": "LAI.html#Estimation-methods-1",
    "page": "LAI Inversion",
    "title": "Estimation methods",
    "category": "section",
    "text": "This package implements 5 different estimation methods, of which we recommend to use ellips_LUT and ellips_opt:Fixed zenith angle of 1 radian ≈ 57.5ᵒ: zenith57\nMiller\'s method miller\nLang\'s method lang\nEllipsoidal leaf angle distribution with Lookup Table estimation ellips_LUT\nEllipsoidal leaf angle distribution with optimization estimation ellips_optAll methods take as argument a PolarImage and a threshold, eg ellips_LUT(polarimg, thresh)."
},

{
    "location": "LAI.html#Zenith-57-1",
    "page": "LAI Inversion",
    "title": "Zenith 57",
    "category": "section",
    "text": "At a viewing angle theta_V of 1 rad ≈ 57.5^o the projection function is almost independent of the leaf angle distribution. For example, with the ellipsoidal distribution we getALIA_to_x(ALIA) = (ALIA/9.65).^-0.6061 - 3\nG(θᵥ, χ) = cos(θᵥ) * sqrt(χ^2 + tan(θᵥ)^2) / (χ+1.702*(χ+1.12)^-0.708)\nusing Winston\np = plot()\nfor θᵥ in linspace(0, pi/2-0.01, 20)        \n    oplot(θL -> G(θᵥ, ALIA_to_x(θL)), 0.1, pi/2, \"--\")     \nend\nθᵥ57 = 1 #in radian\noplot(θL -> G(θᵥ57, ALIA_to_x(θL)),  0.1, pi/2,\"-r\")\ntitle(\"different projection area \\G for a range of view angles\")\nxlabel(\"average leaf inclination angle \\\\theta_L \")(Image: Projection function for different average leaf angles)This method is used for the initial starting point of the ellips_opt method. LAI57 = zenith57(polimg, thresh)"
},

{
    "location": "LAI.html#Miller\'s-method-1",
    "page": "LAI Inversion",
    "title": "Miller\'s method",
    "category": "section",
    "text": "Assuming a constant leaf angle, by integration the effect of G(theta) disappears (Miller 1967):L = 2 int_0^pi2 -ln(P_0(theta))costheta sintheta mathrmd thetaFor this method you need the entire viewing angle up to π/2, which might prove difficult for larger zenith angles (Weiss et al. 2004[1]). Furthermore, the integration over the discrete polar distances of each pixel requires an ambiguous choice of grouping consecutive rings."
},

{
    "location": "LAI.html#Lang\'s-method-1",
    "page": "LAI Inversion",
    "title": "Lang\'s method",
    "category": "section",
    "text": "Lang 1986 approximated -ln(P_0(theta))costheta = K(theta) approx a + btheta around theta_V of 1 rad and obtained L = 2(a+b)We follow Weiss et al. 2004[1] to regress between 25^o and 65^o."
},

{
    "location": "LAI.html#Ellipsoidal-ALIA-estimation-1",
    "page": "LAI Inversion",
    "title": "Ellipsoidal ALIA estimation",
    "category": "section",
    "text": "We can also estimate the parameter chi of the ellipsoidal leaf distribution together with the LAI. We follow Thimonier et al. 2010[2] and use the contact frequency K, with K(theta_V) = G(theta_V chi)L_e = -lnT(theta_v) costheta_Vas the fitting observable. We found more variance over the view zenith range with K (plot below) compared to using ln T(theta_V) (as Norman & Campbell 1989) or T(theta_V) = exp(-G(theta_V chi)L_ecostheta_V) (as Weiss et al. 2004[1]).```julia     ALIA_to_x(ALIA) = (ALIA/9.65).^-0.6061 - 3     G(θᵥ, χ) = cos(θᵥ) * sqrt(χ^2 + tan(θᵥ)^2) / (χ+1.702*(χ+1.12)^-0.708)     using Winstonalia = (5:5:85)*π/180\np = plot()\n\nf(i) = θᵥ ->  G(θᵥ, ALIA_to_x(alia[i]))\n# uncomment below for comparison\n#f(i) = θᵥ ->  exp(-G(θᵥ, ALIA_to_x(alia[i])) / cos(θᵥ))\n#f(i) = θᵥ ->  - G(θᵥ, ALIA_to_x(alia[i])) / cos(θᵥ)\n\noplot(f(1), 0, π/2-0.1, \"-r\")\nfor i = 2:length(alia)-1\n    oplot(f(i), 0, π/2-0.1, \"--\")\nend\noplot(f(int(length(alia)/2)), 0, π/2-0.1, \"-g\") #pi/4\noplot(f(length(alia)), 0, π/2-0.1, \"-k\")\nylabel(\"projection function G\")\nxlabel(\"view angle \\\\theta_V\")(Image: Different fitting variables)Both the Lookup Table approach ellips_LUT(polarimg, thresh) from [Weiss et al. 2004][1] and the optimization method ellips_opt(polarimg, thresh)from Thimonier et al. 2010[2] are implemented. We do not weight the different gap fractions per zenith angle as in Thimonier et al. 2010[2], but we use weighted rings with each a similar amount of pixels. We also use more view zenith rings than in the originals papers because digital cameras have much more pixels these day.We find the parameter space to optimize is smooth and can be seen in a heat map with LUT values and 25 closest solutions in red circles:(Image: Ellipsoidal LUT parameter space)[1]: Weiss et al 2004[2]: Thimonier et al 2010"
},

{
    "location": "clumping.html#",
    "page": "Clumping",
    "title": "Clumping",
    "category": "page",
    "text": ""
},

{
    "location": "clumping.html#Clumping-1",
    "page": "Clumping",
    "title": "Clumping",
    "category": "section",
    "text": "The inversion models assumes a random distribution of the foliage that is often not observed in forest canopies. Usually the outcome of the inversion procedure is considered the effective LAI L_e and a clumping correction Omega is applied to obtain the true LAI L = L_e  Omega. Often Omega is considered dependent on the view zenith angle Omega(theta).We provide two standard methods to calculate the clumping for a specific zenith range: the Lang Xiang method and the Chen Cihlar."
},

{
    "location": "clumping.html#Lang-Xiang-1",
    "page": "Clumping",
    "title": "Lang Xiang",
    "category": "section",
    "text": "The Lang Xiang method divides the zenith ring in azimuthal slices and takes the logarithmic average to calculate the correction (overline denotes the mean):Omega_LX = fracln overline Toverlineln TUse langxiang(polimg, thresh, θ1, θ2, nϕ) to obtain the clumping correction for a zenith range between θ1 and θ2 with nϕ the number of azimuthal slices. "
},

{
    "location": "clumping.html#Chen-Cihlar-1",
    "page": "Clumping",
    "title": "Chen Cihlar",
    "category": "section",
    "text": "The Chen Cihlar method was originally proposed not for hemisperical pictures but for direct (one-dimensional) gap fraction measurements over a line under the canopy. We adopt the method to hemisperical pictures by taking a spiral of a single pixel width, similar to a Fermat\'s spiral, as a basis to calculate the gap fraction (inspiration from private communication with P. Schleppi).The principle for the Chen Cihlar method is to remove large gaps from the measured (cumulative) gap size distribution F_m to get a reduced distribution F_mr until it corresponds to the theoretical (cumulative) gap size distribution F from a random canopy.The theoretical cumulative gap size distribution F depends on the projected LAI L_p and the projected characteristic leaf width W_p:F(lambda theta) = big 1+L_p(theta)lambdaW_p(theta) big expBig(-L_p(theta) -L_p(theta)lambdaW_p(theta) Big)where lambda represents the gap size.W_p can be estimated from the probe probabilty P, i.e. \"the probability that a probe of length l falls completely within a sunfleck on the forect floor\" (Leblanc et al. 2005):W_p = - fracln P(0theta)bigpartial ln P(l theta)  partial lbig_l=0and P(l) = 1L_t sum_i H(lambda_i - l)(lambda_i - l)with lambda_i all the measured gaps, H the Heaviside step function and L_t for normalization.L_p(theta) is initially estimated from the gap fraction T(theta) = -ln T(theta) = -ln F(0theta). Then by reducing F_m to F_mr by excluding those gaps too improbably according the theoretical gap size distributin F, a new value L_p is set by L_p = -ln F_mr and this procedure is repeated a few steps until the reduced distribution matches the theoretical one.The clumping correction is then defined asOmega_CC = fracln F_mln F_mrfrac1-F_mr1-F_m"
},

{
    "location": "calibration.html#",
    "page": "Calibration",
    "title": "Calibration",
    "category": "page",
    "text": ""
},

{
    "location": "calibration.html#Calibration-1",
    "page": "Calibration",
    "title": "Calibration",
    "category": "section",
    "text": "The necessary calibration of the camera & lens setup requires two parts: the lens projection function and the lens center. Both are required to map the pixels of the picture to corresponding polar coordinates of the scene.The lens projection function is assumed to have the same shape for every lens type, but different pixel radius for the 180^o viewing circle of each specific lens with camera setup. This means you can simply take the parameters of the generic projection function for your lens type and multiply by the visible radius. The lens center is assumed different for each different lens with camera, even between same lens type and camera type. However, determining the lens center is straightforward.To determine pixel coordinates of calibration images we recommend to use the ImageView julia package or the free software paint.net."
},

{
    "location": "calibration.html#Lens-Center-1",
    "page": "Calibration",
    "title": "Lens Center",
    "category": "section",
    "text": "The easiest way to determine the lens center is the prick a few holes in the lens cap (or in an attached alumium foile on top of the lens) and use  pictures of the holes turned several times. These holes should each lie on a circle with the lens center in the middle. We provide an IJulia notebook with the necessary straightforward calculations here"
},

{
    "location": "calibration.html#Lens-Projection-Function-1",
    "page": "Calibration",
    "title": "Lens Projection Function",
    "category": "section",
    "text": "The lens projection function maps the view angle θ of the scene to the polar distance ρ to the center. The view angle θ lies in the interval [0, π/2]. The polar distance can be assumed to lie in an interval between [0, 1] for a generic lens projection function, and can then be multiplied by the maximum radius pixels of the specific camera setup. The calibration expects the function for the actual pixels, so multiplied by the maximum radius. The maximum radius can be easily determined by sampling a few points on the horizon of a level test picture and measure the distance to the lens center.Also the inverse projection function is required by calibration method and simply maps back from radius to zenith angle.Note that for both projection functions the calibration method you can provide any function between two scalars in the appropriate domain, not just a polynomial function.We follow the setup from the CAN_EYE user guide and again provide an IJulia notebook with the implementation."
},

]}
