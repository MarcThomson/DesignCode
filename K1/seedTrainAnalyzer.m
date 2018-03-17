ttotal=[0];
Ctotal = zeros(41,1);
VCDtotal =[0];
MabTotal =[0];
NMABTotal = [0];
NCTotal = [0];
O2Total = [0];
%RTotal = [0];
for i = 1:length(vesselSize)
t = totalStructure(-i);C=totalStructure(i);
VCDtemp = VCD(t).*C(11,:);
VCDtotal = [VCDtotal,VCDtemp];
ttotal = [ttotal,ttotal(end)+t];
MabTotal = [MabTotal,C(4,:)];
NCTotal = [NCTotal,VCDtemp*vesselSize(i)];
NMABTotal = [NMABTotal,C(4,:)*vesselSize(i)];
%RTotal = [RTotal,C(end,:)];
Ctotal=[Ctotal,C];
O2Total = [O2Total,C(34,:)*vesselSize(i)/1e6];
end

figure(1);clf;plot(ttotal(2:end),VCDtotal(2:end)/2.31)
figure(2);clf;plot(ttotal(2:end),MabTotal(2:end)*126/1000)
figure(3);clf;plot(ttotal(2:end),NCTotal(2:end)/2.31)
figure(4);clf;plot(ttotal(2:end),NMABTotal(2:end)*126/10^9)
figure(5);plot(ttotal(2:end),Ctotal(end,2:end))
figure(6);plot(ttotal(2:end),O2Total(2:end))


%VCDtotal = VCDtotal(2:end)/2.31;
C = Ctotal(:,2:end);
t = ttotal(2:end);
figure;plotTool;
