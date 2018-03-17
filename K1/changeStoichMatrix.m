initialFluxVec;
stoichMatrix;
internalLogicals;
I = find(internal);
Aeq = stoichMatrix(I,:);


I = find(Aeq);

%for nChanged = 1:2
%    for i = 1:length(I)
        