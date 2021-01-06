function newCentre = findCentre(approxCen,IMG,dMIN,dMAX,nTestPoints)
testCentresX = (approxCen.x-nTestPoints):(approxCen.x+nTestPoints);
testCentresY = (approxCen.y-nTestPoints):(approxCen.y+nTestPoints);
nX = length(testCentresX);
nY = length(testCentresY);
meanErrors = zeros(nY,nX);
for x = 1:nX
    for y = 1:nY
        cen.x = testCentresX(x);
        cen.y = testCentresY(y);
        meanErrors(y,x) = calcErrorsForAssumedCen(cen,IMG,dMIN,dMAX);
    end
end
%imtool(meanErrors);
errorsX = sum(meanErrors,1);
errorsY = sum(meanErrors,2)';

px = polyfit(testCentresX,errorsX,2);
py = polyfit(testCentresY,errorsY,2);
px_values = polyval(px,testCentresX);
py_values = polyval(py,testCentresY);
newCentre.x = -px(2)*0.5/px(1);
newCentre.y = -py(2)*0.5/py(1);
% newCentre
%figure, plot(testCentresX,errorsX,testCentresX,px_values);
%figure, plot(testCentresY,errorsY,testCentresY,py_values);

function meanSqDiff = calcErrorsForAssumedCen(centre,img,dMin,dMax)
q1pixels_y = (centre.y - dMin):-1:(centre.y - dMax);
q1pixels_x = (centre.x-dMin):-1:(centre.x-dMax);
q1 = img(q1pixels_y,q1pixels_x);
q2pixels_y = (centre.y-dMin):-1:(centre.y-dMax);
q2pixels_x = (centre.x+dMin):(centre.x+dMax);
q2 = img(q2pixels_y,q2pixels_x);
q3pixels_y = (centre.y+dMin):(centre.y+dMax);
q3pixels_x = (centre.x-dMin):-1:(centre.x-dMax);
q3 = img(q3pixels_y,q3pixels_x);
q4pixels_y = (centre.y+dMin):(centre.y+dMax);
q4pixels_x = (centre.x-dMin):-1:(centre.x-dMax);
q4 = img(q4pixels_y,q4pixels_x);
IMGvalues(:,:,1) = q1;
IMGvalues(:,:,2) = q2;
IMGvalues(:,:,3) = q3;
IMGvalues(:,:,4) = q4;
errors = std(IMGvalues,0,3);
meanSqDiff = mean(errors(:));
end
end