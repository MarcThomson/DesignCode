ttotal=[0];
Ctotal = zeros(41,1);
VCDtotal =[0];
MabTotal =[0];
NMABTotal = [0];
NCTotal = [0];
%RTotal = [0];
for i = 1:8
t = totalStructure(-i);C=totalStructure(i);
VCDtemp = VCD(t).*C(11,:);
VCDtotal = [VCDtotal,VCDtemp];
ttotal = [ttotal,ttotal(end)+t];
MabTotal = [MabTotal,C(4,:)];
NCTotal = [NCTotal,VCDtemp*vesselSize(i)];
NMABTotal = [NMABTotal,C(4,:)*vesselSize(i)];
%RTotal = [RTotal,C(end,:)];
Ctotal=[Ctotal,C];
end

figure(1);clf;plot(ttotal(2:end),VCDtotal(2:end)/2.31)
figure(2);clf;plot(ttotal(2:end),MabTotal(2:end)*126/1000)
figure(3);clf;plot(ttotal(2:end),NCTotal(2:end)/2.31)
figure(4);clf;plot(ttotal(2:end),NMABTotal(2:end)*126/10^9)
figure(5);plot(ttotal(2:end),Ctotal(end,2:end))