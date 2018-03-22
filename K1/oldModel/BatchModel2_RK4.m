clear;close all;
%% Inputs %%
tend = 10;
cellType = 1;
shift = 0;
Perfusion = 0;
h = 0.1;
shiftDay = 11;
%% Load variables from definitions
reversibleLogicals;
internalLogicals;
initialConditions;
load('stoichMatrix.mat')
t = 0:h:tend;

I = find(~internal);
%C0 = C0(1:end-1);
C0 = initialConditions_vec(I);

C = zeros(length(C0),length(t));
R = zeros(1,length(t));
objVec = zeros(length(t),1);
exitVec = zeros(length(t),1);

%R(1) = C0(end);%
R(1) = initialConditions_vec(end);

C(:,1) = C0;
exitVec(1) = 1;
objVec(1) = -12;
v = zeros(34,length(t));
for i = 1:length(t)-1
    if t(i)>shiftDay
        shift = 1;
        load('stoichMatrix31.mat');
    end
    
   [K1, R1,~,~,~] =  instantRatesV2(C(:,i),t(i),Perfusion,cellType, shift, R(i));
   K1 = K1*h; K1(14) = 0;
   [K2, R2,~,~,~] = instantRatesV2(C(:,i)+K1/2, t(i) + h/2 ,Perfusion,cellType, shift, R1);
   K2 = K2*h; K1(14) = 0;
   [K3, R3,~,~,~] = instantRatesV2(C(:,i)+K2/2, t(i) + h/2 ,Perfusion,cellType, shift, R2);
   K3 = K3*h;K1(14) = 0;
   [K4, R4,~,~,~] = instantRatesV2(C(:,i)+K3, t(i) + h ,Perfusion,cellType, shift, R3);
   K4 = K4*h;  K1(14) = 0;   
   
   rates = (K1+2*K2+2*K3+K4)/6;
   rates(end-1) = 0;
   

   Rtest = (R1+2*R2+2*R3+R4)/6;
   C(:,i+1) = C(:,i) + rates;
   [~, R(i+1),v(:,i+1),objVec(i+1),exitVec(i+1)] =...
       instantRatesV2(C(:,i+1), t(i) + h ,Perfusion,cellType, shift, Rtest);
   
   if  t(i)-floor(t(i))==0 && t(i)>=3
       C(8,i+1) = C(8,i+1) + 8; %GLC
       C(6,i+1) = C(6,i+1)+ 0.1; %CC 
       C(3,i+1) = C(3,i+1) + 2.5  ; %ASN 
       C(9,i+1) = C(9,i+1); %GLN 
       C(15,i+1) = C(15,i+1) + 1; %SER 
   end
  
   
   Ineg = find(C(:,i+1)<0);
   C(Ineg,i+1) = 0;
  
   if abs(R(i+1))>10
       R(i+1) = 0;
       break
   end
end


plotToolV2

figure(2);clf;
for i = 1:34
    subplot(7,5,i);plot(t(2:end),v(i,2:end),'LineWidth',2);
    hold on; plot(t(2:end),0*t(2:end),'r-');xlim([0,10])
    title(num2str(i));
    hold off
end

% figure(3);clf;
% I = find(internal);
% Aeq = stoichMatrix(I,:);
% Aeq(19,end) = 0;
% Aeq(18,end) = 0;
% hold on;
% I1 = find(Aeq(19,1:end-1)>0);
% Aeq1 = Aeq;
% Aeq1(19,I1) = 0;
% plot(t(2:end),-Aeq1(19,:)*v(:,2:end),'b-','LineWidth',2)
% 
% I1 = find(Aeq(19,1:end-1)<0);
% Aeq1 = Aeq;
% Aeq1(19,I1) = 0;
% plot(t(2:end),Aeq1(19,:)*v(:,2:end),'m-','LineWidth',2)
% 
% I1 = find(Aeq(18,1:end-1)<0);
% Aeq1 = Aeq;
% Aeq1(18,I1) = 0;
% plot(t(2:end),Aeq1(18,:)*v(:,2:end),'y-','LineWidth',2)
% 
%  I1 = find(Aeq(18,1:end-1)~=0);
%  Aeq1 = Aeq;
%  Aeq1(18,I1) = 0;
%  plot(t(2:end),Aeq(18,:)*v(:,2:end),'r-','LineWidth',2)

% I1 = find(Aeq(18,1:end-1)>0);
% Aeq1 = Aeq;
% Aeq1(18,I1) = 0;
% plot(t(2:end),Aeq(18,:)*v(:,2:end),'r-','LineWidth',2)
