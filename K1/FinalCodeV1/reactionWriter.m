% Reaction writer
load('StoichMatrix.mat')

names = {'AKG',  '$\red{ALA}$',  'ALA',  '$\red{ANTI}$',  'ANTI',  '$\red{ASN}$',...
    'ASN',  '$\red{ASP}$',  'ASP',  'ATP',  '$\red{BIOM}$',  'BIOM',...
    '$\red{C_C}$',  'C_C',  '$\red{CO_2}$',  '$CO_2$',  'CYS',...
    '$FADH_2$',  'G6P',  '$\red{GLC}$',  '$\red{GLN}$',  'GLN',...
    '$\red{GLU}$',  'GLU',  '$\red{GLY}$',  'GLY',  '$\red{LAC}$',...
    'LAC',  'MAL',  '$\blue{NADH}',  'NADH',  '$\red{NH3}$',  'NH3',...
    '$\red{O_2}$',  '$O_2$',  'OXA',  'PYR',  '$\red{SER}$',  'SER'};
rxnVec = {};

for i1 = 1:length(stoichMatrix(1,:))
    
    I_positive = find(stoichMatrix(:,i1) > 0);
    I_negative = find(stoichMatrix(:,i1) < 0);
    
    % initialize reaction string
    rxn= '';

    % add components consumed to the reaction string
    for i = 1:length(I_negative)-1
        rxn = [rxn,num2str(-stoichMatrix(I_negative(i),i1)),' ',names{I_negative(i)},' + '];
    end

    % add reaction arrow
    rxn = [rxn,num2str(-stoichMatrix(I_negative(end),i1)),' ',names{I_negative(end)},' --> '];

    % add components being produced
    for i = 1:length(I_positive)-1
        rxn = [rxn,num2str(stoichMatrix(I_positive(i),i1)),' ',names{I_positive(i)},' + '];
    end
    rxn = [rxn,num2str(stoichMatrix(I_positive(end),i1)),' ',names{I_positive(end)}];
    rxnVec{i1} = rxn;
end
