function IMG = spotIntensities2img(IntensitiesVec,spotPositions,imgSize,spotDiam)

IMG0 = zeros(imgSize);
nspots = length(IntensitiesVec);
for ii = 1:nspots
    IMG0(spotPositions.y(ii), spotPositions.x(ii)) = IntensitiesVec(ii);
end
IMG = convolveCircle(IMG0,spotDiam);
end
