function displayImage_gray(IMG,invert)

figure();
imagesc(IMG);
pbaspect([size(IMG,2),size(IMG,1),1]);
if invert
    colormap(flipud(gray));
else
    colormap(gray)
end
end