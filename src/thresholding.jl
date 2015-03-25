
threshold(gray) = RidlerCalvard(gray)
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
