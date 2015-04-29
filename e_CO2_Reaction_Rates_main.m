clc 
clear

%Loads the xsec data
file = 'Xsec2.xlsx';

[~,Names,~] = xlsread(file,'B3:AK3'); % the Xsec names
[~,identity,~] = xlsread(file,'B4:AK4'); % the subplot identifiers
[~,Reactant,~] = xlsread(file,'B5:AK5'); %the first reactant, electron or ion
[~,Reactant2,~] = xlsread(file,'B6:AK6'); %the second reactant, what the ion or electron is reacting with
[XsecE,~,~] = xlsread(file,'A7:A67'); % the energies for which the X-secs pertain
[Xsec,~,~] = xlsread(file,'B7:AK67'); % the X-secs, correlated with the energies above

SubPlots = {'Diss';'QuasiDiss';'NonDiss'}; % The identifiers for making plots, must be the same as those used for identifiers
[SubPlotSize,~] = size(SubPlots);

Counter = zeros(SubPlotSize,1);

% Constants (Don't change unless you know what you are doing)
Emin = 0;  % eV, integral energy minimum
Emax = 30; % eV, integral energy maximum
deltaE = .5; % eV, integral step size

%Variables
Ylimit = [10^-12,10^-2]; %The ylimit of the final graph
ArCO2ratio = 1; % 1 is 100% CO2
COperc = .05; % percent, percent of neutral gas which is CO, divide by 2 to get the percent O2
me = 9.10938291e-31; % kg, electron mass
mi = 42*1.67262178*10^-27; % kg, ion mass
RxnTime = (10*10^-9)*30000*2; % s, the streamer cascade time

% density and temperature values used to calculate the reaction rates
global eDensity iDensity nDensity CO2Density ArDensity CODensity O2Density
eDensity = 1e12; % cm^-3, electron density
iDensity = eDensity; % cm^-3, ion density
nDensity = 2.4462670257666400000e19; % cm^-3, neutral density
CO2Density = nDensity*ArCO2ratio*(1-COperc);  % cm^-3
ArDensity = nDensity*(1-ArCO2ratio);  % cm^-3
CODensity = nDensity*ArCO2ratio*COperc;  % CO density, assume 5% conversion
O2Density = CODensity/2;
Te = [.033;.05;.1;.25;.5;1;2;3;4;5;6;7;8;9;10]; % eV, electron temperature 

[L,~] = size(Te);
[~,XsecSize] = size(Xsec);
Ti(1:L,1) = .1;

for X = 1:L
    %Calculate the electron and ion energy distribution functions.  The
    %electrons follow a maxwellian distribution in a DBD plasma and ions
    %follow a druyvestyn distribution.
    [eDist] = VelDist(Emin,Emax,deltaE,Te(X),eDensity,1);
    [iDist] = VelDist(Emin,Emax,deltaE,Ti(X),iDensity,2);
    Counter = zeros(3,1);
    for Y = 1:XsecSize
        %Calculate the reaction rates for the various reactions
        Q=0;
        Z=1;
        while Q == 0
            if strcmp(SubPlots{Z},identity{Y}) == 1
                Counter(Z) = Counter(Z)+1;
                [TotRxnRate(X,Counter(Z),Z)] = RxnTime*RxnRate(Xsec, XsecE,deltaE,Y,eDist,iDist,mi,me,Reactant{Y},Reactant2{Y},ArCO2ratio,COperc);
                if X == 1   
                    RxnRateName(Z,Counter(Z)) = Names(Y);
                end
                Q = 1;
            elseif Z >= SubPlotSize
                fprintf('error in Z \n')
                Q = 1;
            else
                Z=Z+1;
            end
        end
    end
end

%%

%Separates the reaction rates into three graphs
for Q = 1:SubPlotSize
    
    [~,datasize] = size(RxnRateName(Q,:));
    EditRxnRate = TotRxnRate(:,:,Q);
    EditRxnRateName = RxnRateName(Q,:);
    
    %An exception in the third plot that finds the total reaction rate
    if Q == 3
       EditRxnRate(:,datasize+1) = EditRxnRate(:,2)+EditRxnRate(:,4)+EditRxnRate(:,6)+EditRxnRate(:,9);
       EditRxnRateName(datasize+1) = cellstr('Tot Ion');
    end
    R=1;
    S=0;
    % removes reaction rates that fall underneath the rate threshold for
    % graphing, usually they are too low to matter.
    while S == 0
        if max(EditRxnRate(:,R)) < Ylimit(1)
            EditRxnRate(:,R) = [];
            EditRxnRateName(R) = [];
            Counter(Q) = Counter(Q)-1;
            [~,datasize] = size(EditRxnRateName);
        else
            R=R+1;
        end
        if R > datasize
            S = 1;
        end
    end
    % Plot the Rxn Rates and output the data to a text file to be read by gnuplot
    Exporttxt( strcat('../e_CO2_gnuplotdata/',SubPlots{Q}, '.txt'), EditRxnRateName,'T_e',Te,EditRxnRate)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    subplot(1,SubPlotSize,Q)
    semilogy(Te,EditRxnRate)
    set(gca,'fontsize', 12)
    title('e-CO_2 Reaction Rates vs. T_e','fontSize', 14)
    ylim(Ylimit)
    xlabel('T_e (eV)','fontSize', 14)
    ylabel('Reaction Rate (Mol/sec cm^3)','fontSize', 14)
    legend(EditRxnRateName,'fontSize', 9, 'Color', 'none','location','best')
    legend boxoff
end

%pause
%Print(gcf, '-r300', '-dpng','C:\Users\Rufus\Dropbox\Papers\2014\Dissertation paper\latex - Frontiers\images\e_CO2_RxnRates.png') 
%%

% set(0,'DefaultAxesLineStyleOrder','-|--|:|-.','DefaultAxesColorOrder',[1 0 0;0 0 1;1 0 1;0 0 0;0 1 0;0 1 1])
% % Plot the reaction rates versus the electron temperature
% subplot(1,2,1)
% 
% 
% 
% subplot(1,2,2)
% semilogy(Te,NonDiss(:,1:NonDissNameSize-1))
% title('e-CO_2 Reaction Rates vs. T_e')
% ylim(Ylimit)
% xlabel('T_e (eV)','fontSize', 12)
% ylabel('Reaction Rate (Mol/sec*cm^3)','fontSize', 12)
% legend(NonDissName(1:NonDissNameSize-1),'fontSize', 10)
%
% subplot(2,4,3)
% semilogy(Te,O2)
% title('e-O2 Reaction Rates vs. T_e')
% ylim(Ylimit)
% xlabel('T_e (eV)','fontSize', 12)
% ylabel('Reaction Rate (Mol/sec*cm^3)','fontSize', 12)

%%
%
%
% % Same as above except plotting the reaction rates versus the electron
% % density
% Te = 1;
% Ti = .033;
% eDensity = [10^10,10^11,10^12,10^13,10^14,10^15];
%
% [~,L] = size(eDensity);
%
% CO2 = zeros(9,L);
% CON2 = zeros(4,L);
% O2 = zeros(7,L);
%
% for X = 1:L
% %Calculate the electron energy distribution
% [eDist] = VelDist(Emin,Emax,deltaE,Te,eDensity(X));
% [iDist] = VelDist(Emin,Emax,deltaE,Ti,eDensity(X));
%
% %Calculate the reaction rates for the various reactions
% CO2(1,X) = RxnRate(Xsec,deltaE,2,eDist,iDensity, me);  % moles/cm^3*s, CO2 Electron Recombination
% CO2(2,X) = RxnRate(Xsec,deltaE,3,eDist,CO2Density, me); % CO2 Dissociative Attachment
% CO2(3,X) = RxnRate(Xsec,deltaE,4,eDist,CO2Density, me); %CO2 e-Impact ionization
% CO2(4,X) = RxnRate(Xsec,deltaE,5,eDist,CO2Density, me); % CO2 Dissociative Ionization into CO+
% CO2(5,X) = RxnRate(Xsec,deltaE,6,eDist,CO2Density, me); % CO2 Dissociative ionization into O+
% CO2(6,X) = RxnRate(Xsec,deltaE,7,eDist,CO2Density, me); % CO2 e-impact Dissociation
% CO2(7,X) = RxnRate(Xsec,deltaE,8,eDist,CO2Density, me); % CO2 Vibrational Excitation (100)
% CO2(8,X) = RxnRate(Xsec,deltaE,9,eDist,CO2Density, me); % CO2 Vibrational Excitation (010)
% CO2(9,X) = RxnRate(Xsec,deltaE,10,eDist,CO2Density, me); % CO2 Vibrational Excitation (001)
% CON2(1,X) = RxnRate(Xsec,deltaE,11,eDist,CODensity, me); % CO e-impact Ionization
% CON2(2,X) = RxnRate(Xsec,deltaE,12,eDist,CODensity, me); % CO e-impact Dissociation
% CON2(3,X) = RxnRate(Xsec,deltaE,13,iDist,CO2Density, mi); % N2+ + CO2 Charge exchange
% CON2(4,X) = RxnRate(Xsec,deltaE,14,iDist,CO2Density, mi); %N+ + CO2 Charge Exchange
% O2(1,X) = RxnRate(Xsec,deltaE,15,eDist,O2Density, me); % O2 vibrational excitation v'=1
% O2(2,X) = RxnRate(Xsec,deltaE,16,eDist,O2Density, me); % O2 vibrational excitation v'=2
% O2(3,X) = RxnRate(Xsec,deltaE,17,eDist,O2Density, me); % O2 vibrational excitation v'=3
% O2(4,X) = RxnRate(Xsec,deltaE,18,eDist,O2Density, me); % O2 vibrational excitation excitation under 1 eV
% O2(5,X) = RxnRate(Xsec,deltaE,19,eDist,O2Density, me); % O2 e-impact Dissociation
% O2(6,X) = RxnRate(Xsec,deltaE,20,eDist,O2Density, me); % O2 e-impact Ionization
% O2(7,X) = RxnRate(Xsec,deltaE,21,eDist,O2Density, me); % O2 Dissociative Attachment
% end
% 
% % Plot the reaction rates versus the electron temperature
% subplot(2,4,5)
% loglog(eDensity,CO2)
% title('e-CO_2 Reaction Rates vs. N_e')
% ylim(Ylimit)
% xlabel('N_e (cm^3)','fontSize', 12)
% ylabel('Reaction Rate (Mol/sec*cm^3)','fontSize', 12)
% legend(Names(1:9),'fontSize', 10)
% 
% subplot(2,4,6)
% loglog(eDensity,CON2)
% ylim(Ylimit)
% title('e-CO and Charge Exchange Reaction Rates vs. N_e')
% xlabel('N_e (cm^3)','fontSize', 12)
% legend(Names(10:13),'fontSize', 10)
% 
% subplot(2,4,7)
% loglog(eDensity,O2)
% title('e-O2 Reaction Rates vs. N_e')
% ylim(Ylimit)
% xlabel('N_e (cm^3)','fontSize', 12)
% legend(Names(14:20),'fontSize', 10)