function [componentRates, R, v, obj, exit] = instantRatesV2(C,t,Perfusion,cellType,shift,R)
%% Extra inputs
if Perfusion
    VCD = 1;
    carryingCapacity = 90;
else
    VCD = min([1.07 - (0.80/((1 + exp(-0.35*(t - 2.50)))^(1/0.25))),1]);
end

%% Call other parameters
parameterDefinitions;
reversibleLogicals;
internalLogicals;
if shift == 0
    load('stoichMatrix.mat')
else
    load('stoichMatrix31.mat')
end

numRxns = length(stoichMatrix(1,:));
numComponents = length(stoichMatrix(:,1));


%%Modify concentrations
I_external = find(internal == 0);
C_new = zeros(numComponents,1);
C_new(I_external) = C;

%% Optimization Initialization
tol = 1e-10;
options = optimoptions('quadprog','OptimalityTolerance',tol,...
    'MaxIterations',30000,'Algorithm','interior-point-convex',...
    'StepTolerance',1e-12,'Display','off');

Aeq = stoichMatrix;
 for j=flip(1:length(I_external))
     Aeq=[Aeq(1:I_external(j)-1,:);Aeq(I_external(j)+1:end,:)];
 end
beq = zeros(length(Aeq(:,1)),1);

A = zeros(numRxns, numRxns);
for j=flip(1:numRxns)
    if reversible(j)==1
        A=[A(1:j-1,:);A(j+1:end,:)];
    else
        A(j,j) = -1;
    end
end

%make sure it concentrations don't go negative
% I_neg = find(C_new(I_external)<1E-5 );
% I_neg = I_external(I_neg);
% 
% [I_halted_neg,J_halted_neg] = find(stoichMatrix(I_neg,:)<0);
% [I_halted_pos,J_halted_pos] = find(stoichMatrix(I_neg,:)>0);
% 
% A = [A; zeros(length(J_halted_neg), numRxns)];
% for j = 1:length(J_halted_neg)
%     A(end-j+1,J_halted_neg(j)) = 1;
% end
% A = [A; zeros(length(J_halted_pos), numRxns)];
% for j = 1:length(J_halted_pos)
%     A(end-j+1,J_halted_pos(j)) = -1;
% end



b = zeros(length(A(:,1)),1);



y = [C_new; shift;  R; cellType];
inputs = mat2cell(y,ones(1,length(y)),1);
x0 = rateLaws(inputs{:});

if Perfusion
    x0(16) = x0(16)*(1 - C(5)*VCD/2.31  / carryingCapacity );
end



I = find(x0);
H = zeros(numRxns,numRxns);
f = zeros(numRxns,1);    
for j = 1:length(I)
       H(I(j),I(j)) = 2/x0(I(j))^2;
       f(I(j)) = -2/x0(I(j)); 
end


[v, obj, exit] = quadprog(H,f,A,b,Aeq,beq,[],[],[],options);

R = (2*v(1)+0.64*v(16))/v(34);
f = @(x)(tanh(1.1*(x-1.5))+1)/2.*x + (tanh(-1.1*(x-1.5))+1)/2.*0.5  ;
%R = f(R);

componentRates = stoichMatrix* v  * (C_new(11)*VCD/2.31)/1000; 

componentRates = componentRates(I_external);
    
    
end