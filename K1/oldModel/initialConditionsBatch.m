C0=[3.0718
    3.8158
   13.1635
    2.8307
   24.1649
    1.1379
   38.9450
   40.9824
    0.6124
    1.9150
    4.6091
   32.9880
    3.2640
    0.1666
    6.7098]*400/2000;


I = find(~internal);
initialConditionsSeedTrain;
C0 = C0 + initialConditions_vec(I);

C0 = [C0;1.6744];

