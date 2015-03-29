
threshold(gray) = RidlerCalvardfast(gray)
threshold(polim::PolarImage) = threshold(pixels(polim))

function RidlerCalvard(gray)
    max_iter = 100
    tol = 1e-6

    thresh = mean(gray)
    count = 1
    while count < max_iter
        whites = gray .> thresh
        high = mean(gray[whites])
        low = mean(gray[~whites])
        thresh_old = copy(thresh)
        thresh = (high + low)/2
        abs(thresh - thresh_old) < tol && break
        count += 1
    end
    return(thresh)
end

#four times faster
function RidlerCalvardfast(gray)
    max_iter = 100
    tol = 1e-6

    thresh = mean(gray)
    count = 1

    # single pass over input for both high and low mean without temporary array allocation
    function highlowmean(gray, thresh)
        lowmean = highmean = zero(thresh)
        lowcount = highcount = 0
        @inbounds for i = 1:length(gray)
            grayi = gray[i]
            if gray[i] > thresh
                lowmean += grayi
                lowcount += 1
            else
                highmean += grayi
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