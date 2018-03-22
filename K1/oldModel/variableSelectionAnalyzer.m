initialFlux;
reversibleLogicals;
N = 1 ;

features = totalStructure(-N);
fluxes = totalStructure(N);

v_new = initialFluxVec;
v_new(features) = v_new(features)+fluxes';

factor = v_new./initialFluxVec;

factor(find(factor~=1))


A = find(~reversible & factor<0);
if length(A)>0
    fprintf('Failure \n');
end