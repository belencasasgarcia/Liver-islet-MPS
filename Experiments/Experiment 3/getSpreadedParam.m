function pVect = getSpreadedParam(pIn,varargin)
% Collects the min and max value of each parameter and retrieves the
% corresponding parameter vector.

pSize = size(pIn);

if nargin > 1
    useP = varargin{1};
else
    useP = 1:1:pSize(2);
end

nStep = 1;
pVect = zeros([numel(useP)*2 pSize(2)]);

for i = 1 : numel(useP)
    [~, minI] = min(pIn(:,useP(i)));
    [~, maxI] = max(pIn(:,useP(i)));
    pVect(nStep,1:pSize(2)) = pIn(minI,:);
    pVect(nStep+1,1:pSize(2)) = pIn(maxI,:);
    nStep = nStep + 2;
end
