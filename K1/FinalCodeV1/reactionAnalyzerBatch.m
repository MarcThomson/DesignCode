% reactionScraper for Batch Reaction 
% Prints out the net reaction for the batch case and computes the amount of
% MAB Produced


% Molar masses of external components in g/mol
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


% names of the external components
names ={
'ALA','ANTI','ASN','ASP','BIOM','C_C','CO2','GLC','GLN','GLU','GLY',...
'LAC','NH3','O2','SER'
};

% index of the component that the reaction will be normalized to (i.e.,
% this stoich coefficienct will be 1
I_unity = 8;

% compute the difference between the mass input and output for each
% component in mg
if exist('totalFeed')
    massChangeIn = (sum(totalFeed')' + C(:,1)*Volume(1)).*MM; %mg
else
    Volume = vesselSize(vessel)/1e6;
    massChangeIn = (C(:,1)*Volume(1)).*MM;
end
massChangeOut = (C(:,end)*Volume(end)).*MM;  %mg
massChange= massChangeOut - massChangeIn;

% normalize the mass change to the unit index
massChangeScaled = massChange/abs(massChange(I_unity ));

% find the components net produced and consumed
I_positive = find(massChangeScaled > 0);
I_negative = find(massChangeScaled < 0);

% initialize reaction string
rxn= '';

% add components consumed to the reaction string
for i1 = 1:length(I_negative)-1
    rxn = [rxn,num2str(-massChangeScaled(I_negative(i1))),' ',names{I_negative(i1)},' + '];
end

% add reaction arrow
rxn = [rxn,num2str(-massChangeScaled(I_negative(end))),' ',names{I_negative(end)},' --> '];

% add components being produced
for i1 = 1:length(I_positive)-1
    rxn = [rxn,num2str(massChangeScaled(I_positive(i1))),' ',names{I_positive(i1)},' + '];
end

% assume all leftover mass is water and add it to the reaction
rxn = [rxn,num2str(massChangeScaled(I_positive(end))),' ',names{I_positive(end)},' + '];
water = sum(massChangeScaled);
rxn = [rxn,num2str(water),' H2O'];

% find extent with respect to the 
extent = -massChange(I_unity)./massChangeIn(I_unity);

%Amount of MAB produced in kg
MAB_Produced = Volume(end)*C(2,end)*126/1000/1000*1000;
