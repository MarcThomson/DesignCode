%reactionScraper
%

MM = [
89  ;	%	ALA   =  1
126	;	%	ANTI  =  2
132 ;	%	ASN   =  3
133 ;	%	ASP   =  4
152 ;   %	BIOM  =  5
240 ;	%	C_C   =  6
44	;	%	CO2   =  7
180	;	%	GLC   =  8
146 ;	%	GLN   =  9
147 ;	%	GLU   =  10
75  ;	%	GLY   =  11
90	;	%	LAC   =  12
17	;	%	NH3   =  13
32	;	%	O2    =  14
105 ;	%	SER   =  15
   ]; %g/mol


names ={
'ALA','ANTI','ASN','ASP','BIOM','C_C','CO2','GLC','GLN','GLU','GLY',...
'LAC','NH3','O2','SER'
};
I_unity = 8;

massChange = (C(:,end)*Volume(end) - sum(totalFeed')'- C(:,1)*Volume(1)).*MM;  %mg

massChange = massChange/abs(massChange(I_unity ));

I_positive = find(massChange > 0);
I_negative = find(massChange < 0);

rxn= '';

for i = 1:length(I_negative)-1
    rxn = [rxn,num2str(-massChange(I_negative(i))),' ',names{I_negative(i)},' + '];
end
rxn = [rxn,num2str(-massChange(I_negative(end))),' ',names{I_negative(end)},' --> '];
    
for i = 1:length(I_positive)-1
    rxn = [rxn,num2str(massChange(I_positive(i))),' ',names{I_positive(i)},' + '];
end
rxn = [rxn,num2str(massChange(I_positive(end))),' ',names{I_positive(end)},' + '];

water = sum(massChange);
rxn = [rxn,num2str(water),' H2O'];


extent = -(C(I_unity,end)-C(I_unity,1))/C(I_unity,1);

MAB_Producted = Volume(end)*C(2,end)*126/1000/1000*1000;
