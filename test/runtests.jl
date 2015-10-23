using LeafAreaIndex
using Base.Test


## CALIBRATION ##

testimg1 = reshape(collect(1:25), 5, 5)
Rmax1 = 2
test1cal = calibrate(5, 5, 3, 3, θ -> 2 * Rmax1 * θ / π, ρ -> ρ * π/2 / Rmax1)

@test_approx_eq test1cal.fθρ(0) 0
@test_approx_eq test1cal.fρθ(0) 0
@test_approx_eq test1cal.fθρ(pi/2) Rmax1
@test_approx_eq test1cal.fρθ(2) π/2
@test test1cal.size1 == size(testimg1, 1)
@test test1cal.size2 == size(testimg1, 2)

@test test1cal.sort_ind == testimg1[test1cal.sort_ind]
@test testimg1[test1cal.spiral_ind] == [13,9,14,19,18,17,12,7,8,15,23,11,3]
@test test1cal.ρ²unique == [0,1,2,4]
@test test1cal.ρ²unique_ind == UInt32[1,2,6,10]

@test test1cal.ρ²sort == UInt32[0,1,1,1,1,2,2,2,2,4,4,4,4]
@test_approx_eq test1cal.ϕsort [0, -π/2:π/2:π, -3π/4:π/2:3π/4, -π/2:π/2:π]


## THRESHOLDING ##

testimg2 = 0.5 * ones(FixedPointNumbers.Ufixed16, 5, 5)
testimg2[2:4, 2:4] = 0.7
testimg2[3, 3] = 1
@test minimum_threshold(testimg2) > 0.5

# CircQueue
let
	a = [1, 2, 3]
	c = CircQueue(a)	
	@test pushshift!(c, 5) == 1
	@test pushshift!(c, 6) == 2
	@test pushshift!(c, 7) == 3
	@test pushshift!(c, 8) == 5
end
