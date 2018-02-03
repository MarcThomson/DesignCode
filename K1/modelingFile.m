%Modeling File
%% input
clear all;close all;
cellType = 1;
dt = 0.25; %days;
tend = 5; %days;
%%
%load variables
parameterDefinitions;
reversibleLogicals;
initialConditions;
internalLogicals;
load('stoichMatrix.mat')

t = 0:dt:tend;
VCD = @(t1)(1-1/(1+exp(-0.5*(t1-9))));

%tol = 1e-15;
%options = optimoptions('fmincon','ConstraintTolerance',tol,'StepTolerance',tol,'OptimalityTolerance',tol);
options = optimoptions('fmincon','MaxFunctionEvaluations',30000);

v = zeros(length(stoichMatrix(1,:)),length(t));
C = zeros(length(stoichMatrix(:,1))+2,length(t)); %includes R and shift
C(:,1) = initialConditions_vec;
 for i=1:length(t)-1
     if i==1
         y = [initialConditions_vec; cellType];
     else
         y = [C(:,i-1); cellType];
     end

     inputs = mat2cell(y,ones(1,length(y)),1);
     rates = rateLaws(inputs {:});
     x0 = rates;

     %linear ops
     I = find(rates~=0);
     Aeq = stoichMatrix;
     beq = zeros(length(Aeq(:,1)),1);
     I_internal=find(internal);
     for j=1:length(I_internal)
         Aeq(I_internal(j),:) = 0;
     end

     A = zeros(length(x0),length(x0));
     for j=1:length(x0)
         A(j,i) = reversible(j)-1;
     end
     b = reversible; 

     f_obj=@(x) sum( ((x(I)-x0(I))./(x0(I))).^2);
     v(:,i+1) = fmincon(f_obj,x0,A,b,Aeq,beq,[],[],[],options);
     R = (2*v(1,i+1)+0.64*v(16,i+1))/v(34,i+1);

     %update new values
     C(1:end-2,i+1) =  C(1:end-2,i) +  (stoichMatrix*v(:,i+1))*dt* VCD(t(i))/1000;
     C(end-1,i+1) = 0;
     C(end, i+1) = R;   
 end