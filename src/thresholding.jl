
threshold(gray) = RidlerCalvard(gray)
threshold(polim::PolarImage) = threshold(pixels(polim))

function RidlerCalvard(gray)
    max_iter = 100
    tol = 1e-6

    @show thresh = mean(gray)
    count = 1

    # Single pass over input for both high and low mean is
    # almost 5 times faster if no temporary array allocation.
    function highlowmean(gray, thresh)
        lowmean = highmean = zero(thresh)
        lowcount = highcount = 0
        for pixel in gray            
            if pixel < thresh
                lowmean += pixel
                lowcount += 1
            else
                highmean += pixel
                highcount += 1
            end
        end
        highmean/highcount, lowmean/lowcount
    end

    while count < max_iter        
        high, low  = highlowmean(gray, thresh)
        thresh_old = thresh
        thresh = (high + low)/2
        abs(thresh - thresh_old) < tol && break
        count += 1
    end
    return(thresh)
end
