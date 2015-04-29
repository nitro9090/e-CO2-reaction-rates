function [Dist] = VelDist(Emin,Emax,deltaE,T,density,Distribution)
%Calculating the velocity distribution
%Calculating the standard
k = 1.6022*10^-19; % J/eV

Xmax = ceil((Emax-Emin)/deltaE)+1;

StdDist = 0;
Velocity = 0;
for X = 1:Xmax
    E = Emin + deltaE*(X-1);
    if Distribution == 1
        Dist1(X) = exp(-1*(E/T)^2); % Druyvesteyn Distribution
    elseif Distribution == 2
        Dist1(X) = exp(-1*E/T); % Maxwellian Distribution
    else
        fprintf('Distribution fail')
    end
    if X > 1
        StdDist(X-1,1) = (Dist1(X-1)+Dist1(X))*deltaE/2;
        Dist(X-1,1) = E-deltaE/2;
    end
end

SumDist = sum(StdDist);
Dist(:,2) = StdDist / SumDist * density;
end