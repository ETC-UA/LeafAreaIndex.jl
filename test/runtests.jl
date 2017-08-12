using LeafAreaIndex
using Base.Test
import FixedPointNumbers

## CALIBRATION ##

testimg1 = reshape(collect(1:25), 5, 5)
Rmax1 = 2
test1cal = CameraLens(5, 5, 3, 3, θ -> 2 * Rmax1 * θ / π, ρ -> ρ * π/2 / Rmax1)

@test test1cal.fθρ(0) ≈ 0
@test test1cal.fρθ(0) ≈ 0
@test test1cal.fθρ(pi/2) ≈ Rmax1
@test test1cal.fρθ(2) ≈ π/2
@test test1cal.height == size(testimg1, 1)
@test test1cal.width == size(testimg1, 2)

@test test1cal.sort_ind == testimg1[test1cal.sort_ind]
@test testimg1[test1cal.spiral_ind] == [13,9,14,19,18,17,12,7,8,15,23,11,3]

@test test1cal.ρ²sort == UInt32[0,1,1,1,1,2,2,2,2,4,4,4,4]
#@test test1cal.ϕsort ≈ collect([0, -π/2:π/2:π, -3π/4:π/2:3π/4, -π/2:π/2:π])

## THRESHOLDING ##

testimg2 = 0.5 * ones(FixedPointNumbers.Normed{UInt16, 16}, 5, 5)
testimg2[2:4, 2:4] = 0.7
testimg2[3, 3] = 1
@test threshold(testimg2, MinimumThreshold()) > 0.5

# circqueue
let
	a = [1; 2; 3]
	c = LeafAreaIndex.circqueue(a)	
	@test LeafAreaIndex.pushshift!(c, 5) == 1
	@test LeafAreaIndex.pushshift!(c, 6) == 2
	@test LeafAreaIndex.pushshift!(c, 7) == 3
	@test LeafAreaIndex.pushshift!(c, 8) == 5
end
