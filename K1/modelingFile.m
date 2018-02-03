%Modeling File
%% input

cellType = 1;

%%

%load variables
parameterDefinitions;
reversibleLogicals;
initialConditions;
internalLogicals;
load('stoichMatrix.mat')

%get initial fluxes;
y = [initialConditions_vec; cellType];
inputs = mat2cell(y,ones(1,length(y)),1);
rates = rateLaws(inputs {:});

x0 = rates;

I = find(rates~=0);

f_obj=@(x) sum( ((x(I)-x0(I))./(x0(I))).^2);

Aeq = stoichMatrix;
beq = zeros(length(Aeq(:,1)),1);

A = zeros(length(x0),length(x0));
for i=1:length(x0)
    A(i,i) = reversible(i)-1;
end
b = reversible; 

tol = 1e-15;
options = optimoptions('fmincon','ConstraintTolerance',tol,'StepTolerance',tol,'OptimalityTolerance',tol);
v = fmincon(f_obj,x0,A,b,Aeq,beq,[],[],[],options);
2+2
% for i=1:1
%     
%     
% 
%     
% end