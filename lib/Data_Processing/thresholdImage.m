function IMG2 = thresholdImage(IMG,threshold)

IMG2 = IMG;
IMG2(find(IMG < threshold)) = 0;
end