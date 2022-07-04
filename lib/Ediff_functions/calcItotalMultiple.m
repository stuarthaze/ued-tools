function [Itotal,Isingle,Imultiple] = calcItotalMultiple(ImultiSeparate,xThick)
% (c) S.Hayes - June 2020
% Calculates total radial scattering with single and multiple scattering
% Input:
%   ImultiSeparate - Size: [nEvents,ndata]
%                  - Output from calcRadialMultipleScattering()
%   xThick         - Sample thickness in terms of elastic mean-free path

nEvents = size(ImultiSeparate,1);
if nEvents < xThick*3
    disp('Warning: Too few multiple scattering events for sample thickness')
end

ImultiWeighted = ImultiSeparate;
for ii = 1:nEvents
    ImultiWeighted(ii,:) = ImultiSeparate(ii,:)*poisspdf(ii,xThick);
end
ImultiWeighted = ImultiWeighted/poisspdf(1,xThick);
Itotal = sum(ImultiWeighted,1);
Isingle = ImultiWeighted(1,:);
Imultiple = Itotal - Isingle;
end
