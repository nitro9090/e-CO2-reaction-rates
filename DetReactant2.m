function [ Reactant2Density] = DetReactant2( Reactant )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global iDensity nDensity CO2Density ArDensity CODensity O2Density

if strcmp(Reactant,'iDensity') == 1
    Reactant2Density = iDensity;
elseif strcmp(Reactant,'CO2Density') == 1
    Reactant2Density = CO2Density; 
elseif strcmp(Reactant,'O2Density') == 1
    Reactant2Density = O2Density;
elseif strcmp(Reactant,'CODensity') == 1
    Reactant2Density = CODensity;
elseif strcmp(Reactant,'ArDensity') == 1
    Reactant2Density = ArDensity;
elseif strcmp(Reactant,'nDensity') == 1
    Reactant2Density = nDensity;
elseif strcmp(Reactant,'CO2IonDensity') == 1
    Reactant2Density = CO2Density/nDensity*iDensity;
elseif strcmp(Reactant,'ArIonDensity') == 1
    Reactant2Density = ArDensity/nDensity*iDensity;
elseif strcmp(Reactant,'COIonDensity') == 1
    Reactant2Density = CODensity/nDensity*iDensity;
elseif strcmp(Reactant,'O2IonDensity') == 1
    Reactant2Density = O2Density/nDensity*iDensity;
else
    fprintf('error with %s\n',Reactant)
    Reactant2Density = 0;
end

end

