function R = partials2cartesians(XYZ,AXES)
TRFRM(1,:) = AXES.a;
TRFRM(2,:) = AXES.b;
TRFRM(3,:) = AXES.c;
R = XYZ*TRFRM;
end