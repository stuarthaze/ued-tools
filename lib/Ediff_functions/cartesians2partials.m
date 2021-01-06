function uvw = cartesians2partials(xyz, AXES)
TRFRM(1,:) = AXES.a;
TRFRM(2,:) = AXES.b;
TRFRM(3,:) = AXES.c;
INVRS = inv(TRFRM);
uvw = xyz*INVRS;
end

