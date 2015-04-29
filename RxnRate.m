function [TotRxnRate] = RxnRate(Xsec,XsecE,deltaE,Pos,eDist,iDist,mi,me,Reactant,Reactant2,ArCO2ratio,COperc)
%Calculates the reaction rates based on the input values

k = 1.6022*10^-19; % J/eV
Fail = 0;

if strcmp(Reactant,'electron')
    Dist = eDist;
    mass = me;
elseif strcmp(Reactant,'ion')
    Dist = iDist;
    mass = mi;
elseif strcmp(Reactant,'ArIon')
    Dist = iDist*(1-ArCO2ratio)*(1-COperc);
    mass = mi;
elseif strcmp(Reactant,'CO2Ion')
    Dist = iDist*(ArCO2ratio)*(1-COperc);
    mass = mi;
else
    fprintf('Exterminate! No density is setup for %s, check the Rxn Rate function \n',Reactant)
    Fail = 1;
end

if Fail == 0
    [Reactant2Density] = DetReactant2(Reactant2);
    
    [Ysize,~] = size(Dist);
    
    for A = 1:Ysize
        X = 0;
        J = 0;
        while X == 0
            J=J+1;
            % if the distribution energy values fall between two tabulated X-sec
            % energies, then it linearly scales the X-sec values to the
            % appropriate energy, otherwise it uses the appropriate X-sec.
            if Dist(A,1) <= XsecE(J)
                X = 1;
                Q = 0;
                MinCt = 0;
                while Q == 0
                    MinCt = MinCt+1;
                    if MinCt >= J
                        MinXsec = 0;
                        Q = 1;
                    elseif Xsec(J-MinCt,Pos) < 0
                    else
                        MinXsec = Xsec(J-MinCt,Pos);
                        Q = 1;
                    end
                end
                
                R = 0;
                MaxCt = -1;
                while R == 0
                    MaxCt = MaxCt+1;
                    if MaxCt >= Ysize-J
                        MaxXsec = 0;
                        R = 1;
                    elseif Xsec(J+MaxCt,Pos) < 0
                    else
                        MaxXsec = Xsec(J+MaxCt,Pos);
                        R = 1;
                    end
                end
                if MaxXsec == MinXsec
                    TotXsec(A) = MaxXsec;
                else
                    TotXsec(A) = MinXsec + (MaxXsec-MinXsec)*(MinCt-1+rem(Dist(A,1),deltaE)/deltaE)/(MaxCt+MinCt);
                end
            end
        end
        Velocity(A) = (2*Dist(A,1)*k/mass)^.5*100; % cm/s, calculating the particle velocity
        RxnRates(A) = Dist(A,2)*TotXsec(A)*Velocity(A)*Reactant2Density; % calculating the rxn rates
    end
    TotRxnRate = sum(RxnRates) / (6.022*10^23); % calculates the total rxn rate for the
else
    TotRxnRate = 0;
end
end