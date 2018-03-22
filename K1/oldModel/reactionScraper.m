%reactionScraper
%
batch = 0; %1 = perfusion

MM = [
89	;	%	ALA   =  1
126	;	%	ANTI  =  2
132	;	%	ASN   =  3
133	;	%	ASP   =  4
152 ;   %	BIOM  =  5
240	;	%	C_C   =  6
44	;	%	CO2   =  7
180	;	%	GLC   =  8
146 ;	%	GLN   =  9
147	;	%	GLU   =  10
75	;	%	GLY   =  11
90	;	%	LAC   =  12
17	;	%	NH3   =  13
32	;	%	O2    =  14
105	;	%	SER   =  15
   ]; %g/mol

if batch == 1
    I_external = [
        2, 4, 6, 8, 11, 13, 15, 20, 21, 23, 25, 27, 32, 34, 38
    ]; 
else
    I_external = (1:15);
end

names ={
'ALA','ANTI','ASN','ASP','BIOM','C_C','CO2','GLC','GLN','GLU','GLY',...
'LAC','NH3','O2','SER'
};

massChange = (C(end,I_external)-C(1,I_external)).*MM';
massChange = massChange/massChange(2);
