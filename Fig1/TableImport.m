clear;
clc;
T1=readtable('D:/PlasmidChromosome/NCBIPlasmids.xlsx','VariableNamingRule','preserve');
T2=readtable('prokaryotes.csv','VariableNamingRule','preserve');
save('DataSets.mat','T1','T2');