function R = partials2cartesians_reciprocal(HKL,AXES)
TRFRM(1,:) = AXES.a_star;
TRFRM(2,:) = AXES.b_star;
TRFRM(3,:) = AXES.c_star;
R = HKL*TRFRM;
end