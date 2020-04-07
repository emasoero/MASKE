clear all
close all



% First hypothesis of reactions
% C3S + 5H2O = 3Ca + H4SiO4 + 6OH      with lnK = -50.7    from Nicoleau 2013
%
% CSHI = (Ca + H3SiO4 + OH + 4H2O) with lnK = -22.62     from Bullard 2015
%       + H3SiO4 + H2O = H4SiO4 + OH with lnK = -9.6   from Bullard 2015
% which gives
% CSHI = Ca + H4SiO4 +2OH +3H2O   with lnK = -32.22
%
% CSHII = 2Ca + H3SiO4 + 3OH  +3H2O   with lnK = -34.54
%       + H3SiO4 + H2O = H4SiO4 + OH with lnK = -9.6   from Bullard 2015
% which gives
% CSHII = 2Ca + H4SiO4 +4OH 2H2O   with lnK = -44.14
%
%
%
% Second hypothesis of reaction, using Bullard 2015 directly
% C3S + 4H2O = 3Ca + H3SiO4 + 5OH      with lnK = -41.1 
%
% CSHI = (Ca + H3SiO4 + OH + 4H2O) with lnK = -22.62     from Bullard 2015
%
% CSHII = 2Ca + H3SiO4 + 3OH  +3H2O   with lnK = -34.54
%



% read experimental solubility data for CSH from Jennings 1986
%fileID = fopen('Jennings86.txt');
%mydata = textscan(fileID,'%f%f','HeaderLines',1);
%fclose(fileID);

% Debye Huckel coefficients
A = 0.51; % mol-1/2 L-1
B = 3.29;  % nm-1 mol-1/2 L-1




% Compute solubility curves for H4SiO4 case

figure
subplot(1,3,1);

% plot Jennings 1986 data first
%semilogy(mydata{1,1},mydata{1,2},'-*')
%hold on

% C3S with H4SiO4
% set target equilibrium constants
Kc3s = exp(-50.7);  %from Nicoleau 2013   -- EDGE DEFECT -- 

%charges
zCa = 2;
zSi = 0;
zOH = -1;

% hydrated radii of molecules
dCa = 0.82;   %nm
dSi = 0.983;   %computed for H4SiO4, not for H3SiO4, but I guess they are similar enough
dOH = 0.416;


% create vector of molar concentrations in solution
vCa = linspace(0.1,25,100)*1e-3;  % this is fixed; the others computed
vSi = linspace(-1,-1,100);
vOH = linspace(-1,-1,100);

% for each Ca concentration in vector
for i=1:length(vCa)
    Ca = vCa(i);
    Si = 0;  % initial assumption of no Si
    OH = 2*Ca; % OH is linked to Ca and Si content to keep charge neutrality
    
    %initial assumption of activity coefficients
    gCa = 1;   gSi = 1;   gOH = 1;
    
    % thus initial activities and their product
    aCa = gCa*Ca;     aSi = gSi*Si;    aOH = gOH*OH;
    Q = (aCa^3)*(aSi)*(aOH^6);  % having set Si to zero, necessariliy leads to Q = 0
    
    %   iterative loop
    tol =  0.01;  % tolerance for convergence
    while abs (Kc3s - Q)/Kc3s > tol
        
        % update concentrations of Si and OH to get Q closer to Kc3s
        Qnew = Q + (Kc3s - Q)/2;
        Si = Qnew/(aCa^3)/(aOH^6)/gSi;
        OH = 2*Ca;   % in general, OH could be function of Si, in which case it needs to be updated here   
        
        % compute activity coefficients g via Debye Huckel
        I = 1/2*(Ca*zCa*zCa + Si*zSi*zSi + OH*zOH*zOH); %ionic strength
        gCa = 10^(-A*zCa*zCa*sqrt(I)/(1+B*dCa*sqrt(I))); %coefficients
        gSi = 10^(-A*zSi*zSi*sqrt(I)/(1+B*dSi*sqrt(I)));
        gOH = 10^(-A*zOH*zOH*sqrt(I)/(1+B*dOH*sqrt(I)));
        
        % compute activities and their product
        aCa = gCa * Ca;
        aSi = gSi * Si;
        aOH = gOH * OH;
        Q = (aCa^3)*(aSi)*(aOH^6);
    end
    vSi(i) = Si;
    vOH(i) = OH;
end


semilogy(vCa*1e3,vSi*1e6,'LineWidth',2)
hold on





% C3S with H4SiO4
% set target equilibrium constants
Kc3s = exp(-63);  %from Nicoleau 2013 -- SCREW DISLOCATION -- 

%charges
zCa = 2;
zSi = 0;
zOH = -1;

% hydrated radii of molecules
dCa = 0.82;   %nm
dSi = 0.983;   %computed for H4SiO4, not for H3SiO4, but I guess they are similar enough
dOH = 0.416;


% create vector of molar concentrations in solution
vCa = linspace(0.1,25,100)*1e-3;  % this is fixed; the others computed
vSi = linspace(-1,-1,100);
vOH = linspace(-1,-1,100);

% for each Ca concentration in vector
for i=1:length(vCa)
    Ca = vCa(i);
    Si = 0;  % initial assumption of no Si
    OH = 2*Ca; % OH is linked to Ca and Si content to keep charge neutrality
    
    %initial assumption of activity coefficients
    gCa = 1;   gSi = 1;   gOH = 1;
    
    % thus initial activities and their product
    aCa = gCa*Ca;     aSi = gSi*Si;    aOH = gOH*OH;
    Q = (aCa^3)*(aSi)*(aOH^6);  % having set Si to zero, necessariliy leads to Q = 0
    
    %   iterative loop
    tol =  0.01;  % tolerance for convergence
    while abs (Kc3s - Q)/Kc3s > tol
        
        % update concentrations of Si and OH to get Q closer to Kc3s
        Qnew = Q + (Kc3s - Q)/2;
        Si = Qnew/(aCa^3)/(aOH^6)/gSi;
        OH = 2*Ca;   % in general, OH could be function of Si, in which case it needs to be updated here   
        
        % compute activity coefficients g via Debye Huckel
        I = 1/2*(Ca*zCa*zCa + Si*zSi*zSi + OH*zOH*zOH); %ionic strength
        gCa = 10^(-A*zCa*zCa*sqrt(I)/(1+B*dCa*sqrt(I))); %coefficients
        gSi = 10^(-A*zSi*zSi*sqrt(I)/(1+B*dSi*sqrt(I)));
        gOH = 10^(-A*zOH*zOH*sqrt(I)/(1+B*dOH*sqrt(I)));
        
        % compute activities and their product
        aCa = gCa * Ca;
        aSi = gSi * Si;
        aOH = gOH * OH;
        Q = (aCa^3)*(aSi)*(aOH^6);
    end
    vSi(i) = Si;
    vOH(i) = OH;
end


semilogy(vCa*1e3,vSi*1e6,'LineWidth',2)
hold on









% CSHI with H4SiO4 --- CSHI = Ca + H4SiO4 +2OH +3H2O
% set target equilibrium constants
KcshI = exp(-32.22);  %from Bullard 2013

% charges
zCa = 2;
zSi = 0;
zOH = -1;

% hydrated radii of molecules
dCa = 0.82;   %nm
dSi = 0.983;   %computed for H4SiO4, not for H3SiO4, but I guess they are similar enough
dOH = 0.416;


% create vector of molar concentrations in solution
vCa = linspace(0.1,25,100)*1e-3;  % this is fixed; the others computed
vSi = linspace(-1,-1,100);
vOH = linspace(-1,-1,100);

% for each Ca concentration in vector
for i=1:length(vCa)
    Ca = vCa(i);
    Si = 0;  % initial assumption of no Si
    OH = 2*Ca; % OH is linked to Ca and Si content to keep charge neutrality
    
    %initial assumption of activity coefficients
    gCa = 1;   gSi = 1;   gOH = 1;
    
    % thus initial activities and their product
    aCa = gCa*Ca;     aSi = gSi*Si;    aOH = gOH*OH;
    Q = (aCa)*(aSi)*(aOH^2);  % having set Si to zero, necessariliy leads to Q = 0
    
    %   iterative loop
    tol =  0.01;  % tolerance for convergence
    while abs (KcshI - Q)/KcshI > tol
        
        % update concentrations of Si and OH to get Q closer to Kc3s
        Qnew = Q + (KcshI - Q)/2;
        Si = Qnew/(aCa)/(aOH^2)/gSi;
        OH = 2*Ca;   % in general, OH could be function of Si, in which case it needs to be updated here   
        
        % compute activity coefficients g via Debye Huckel
        I = 1/2*(Ca*zCa*zCa + Si*zSi*zSi + OH*zOH*zOH); %ionic strength
        gCa = 10^(-A*zCa*zCa*sqrt(I)/(1+B*dCa*sqrt(I))); %coefficients
        gSi = 10^(-A*zSi*zSi*sqrt(I)/(1+B*dSi*sqrt(I)));
        gOH = 10^(-A*zOH*zOH*sqrt(I)/(1+B*dOH*sqrt(I)));
        
        % compute activities and their product
        aCa = gCa * Ca;
        aSi = gSi * Si;
        aOH = gOH * OH;
        Q = (aCa)*(aSi)*(aOH^2);
    end
    vSi(i) = Si;
    vOH(i) = OH;
end


semilogy(vCa*1e3,vSi*1e6)
hold on





% CSHII with H4SiO4 --- CSHII = 2Ca + H4SiO4 +4OH 2H2O
% set target equilibrium constants
KcshII = exp(-44.14);  %from Bullard 2013

% charges
zCa = 2;
zSi = 0;
zOH = -1;

% hydrated radii of molecules
dCa = 0.82;   %nm
dSi = 0.983;   %computed for H4SiO4, not for H3SiO4, but I guess they are similar enough
dOH = 0.416;


% create vector of molar concentrations in solution
vCa = linspace(0.1,25,100)*1e-3;  % this is fixed; the others computed
vSi = linspace(-1,-1,100);
vOH = linspace(-1,-1,100);

% for each Ca concentration in vector
for i=1:length(vCa)
    Ca = vCa(i);
    Si = 0;  % initial assumption of no Si
    OH = 2*Ca; % OH is linked to Ca and Si content to keep charge neutrality
    
    %initial assumption of activity coefficients
    gCa = 1;   gSi = 1;   gOH = 1;
    
    % thus initial activities and their product
    aCa = gCa*Ca;     aSi = gSi*Si;    aOH = gOH*OH;
    Q = (aCa^2)*(aSi)*(aOH^4);  % having set Si to zero, necessariliy leads to Q = 0
    
    %   iterative loop
    tol =  0.01;  % tolerance for convergence
    while abs (KcshII - Q)/KcshII > tol
        
        % update concentrations of Si and OH to get Q closer to Kc3s
        Qnew = Q + (KcshII - Q)/2;
        Si = Qnew/(aCa^2)/(aOH^4)/gSi;
        OH = 2*Ca;   % in general, OH could be function of Si, in which case it needs to be updated here   
        
        % compute activity coefficients g via Debye Huckel
        I = 1/2*(Ca*zCa*zCa + Si*zSi*zSi + OH*zOH*zOH); %ionic strength
        gCa = 10^(-A*zCa*zCa*sqrt(I)/(1+B*dCa*sqrt(I))); %coefficients
        gSi = 10^(-A*zSi*zSi*sqrt(I)/(1+B*dSi*sqrt(I)));
        gOH = 10^(-A*zOH*zOH*sqrt(I)/(1+B*dOH*sqrt(I)));
        
        % compute activities and their product
        aCa = gCa * Ca;
        aSi = gSi * Si;
        aOH = gOH * OH;
        Q = (aCa^2)*(aSi)*(aOH^4);
    end
    vSi(i) = Si;
    vOH(i) = OH;
end


semilogy(vCa*1e3,vSi*1e6)
legend({'Jennings1986','C3S with K = -50.7','C3S with K = -63','CSHI','CSHII'},'FontSize',14)
title('Solubilities assuming only H4SiO4 in solution','FontSize',18)
xlabel('Ca (mmol)', 'FontSize', 18);
ylabel('Si (umol)', 'FontSize', 18);
set(gca,'FontSize',14)
xlim([0 27])
ylim([1e-5 1e10])







% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************



% Compute solubility curves for H3SiO4 case

subplot(1,3,2);

% plot Jennings 1986 data first
%semilogy(mydata{1,1},mydata{1,2},'-*')
%hold on

    

% C3S with H3SiO4     --> C3S + 4H2O = 3Ca + H3SiO4 + 5OH  with lnK = -41.1 
% set target equilibrium constants
Kc3s = exp(-41.1);  %from Nicoleau 2013   -- EDGE DEFECT -- 

%charges
zCa = 2;
zSi = -1;
zOH = -1;

% hydrated radii of molecules
dCa = 0.82;   %nm
dSi = 0.983;   %computed for H4SiO4, not for H3SiO4, but I guess they are similar enough
dOH = 0.416;


% create vector of molar concentrations in solution
vCa = linspace(0.1,25,100)*1e-3;  % this is fixed; the others computed
vSi = linspace(-1,-1,100);
vOH = linspace(-1,-1,100);

% for each Ca concentration in vector
for i=1:length(vCa)
    Ca = vCa(i);
    Si = 0;  % initial assumption of no Si
    OH = 2*Ca-Si; % OH is linked to Ca and Si content to keep charge neutrality
    
    %initial assumption of activity coefficients
    gCa = 1;   gSi = 1;   gOH = 1;
    
    % thus initial activities and their product
    aCa = gCa*Ca;     aSi = gSi*Si;    aOH = gOH*OH;
    Q = (aCa^3)*(aSi)*(aOH^5);  % having set Si to zero, necessariliy leads to Q = 0
    
    %   iterative loop
    tol =  0.01;  % tolerance for convergence
    count = 0;
    while abs (Kc3s - Q)/Kc3s > tol
        
        count = count + 1;
        % update concentrations of Si and OH to get Q closer to Kc3s
        Qnew = Q + (Kc3s - Q)/100;
        if count == 100000
            Si = -1;
            OH = -1;
            break
        end
        Si = Qnew/(aCa^3)/(aOH^5)/gSi;
        OH = 2*Ca-Si;   % in general, OH could be function of Si, in which case it needs to be updated here   
        
        % compute activity coefficients g via Debye Huckel
        I = 1/2*(Ca*zCa*zCa + Si*zSi*zSi + OH*zOH*zOH); %ionic strength
        gCa = 10^(-A*zCa*zCa*sqrt(I)/(1+B*dCa*sqrt(I))); %coefficients
        gSi = 10^(-A*zSi*zSi*sqrt(I)/(1+B*dSi*sqrt(I)));
        gOH = 10^(-A*zOH*zOH*sqrt(I)/(1+B*dOH*sqrt(I)));
        
        % compute activities and their product
        aCa = gCa * Ca;
        aSi = gSi * Si;
        aOH = gOH * OH;
        Q = (aCa^3)*(aSi)*(aOH^5);
    end
    vSi(i) = Si;
    vOH(i) = OH;
end


semilogy(vCa*1e3,vSi*1e6,'LineWidth',2)
hold on





% C3S with H3SiO4     --> C3S + 4H2O = 3Ca + H3SiO4 + 5OH  with lnK = -41.1 
% set target equilibrium constants
Kc3s = exp(-41.1-13);  % approximate estimation by generalising Nicoleau -- SCREW DISLOCATION -- 

%charges
zCa = 2;
zSi = -1;
zOH = -1;

% hydrated radii of molecules
dCa = 0.82;   %nm
dSi = 0.983;   %computed for H4SiO4, not for H3SiO4, but I guess they are similar enough
dOH = 0.416;


% create vector of molar concentrations in solution
vCa = linspace(0.1,25,100)*1e-3;  % this is fixed; the others computed
vSi = linspace(-1,-1,100);
vOH = linspace(-1,-1,100);

% for each Ca concentration in vector
for i=1:length(vCa)
    Ca = vCa(i);
    Si = 0;  % initial assumption of no Si
    OH = 2*Ca-Si; % OH is linked to Ca and Si content to keep charge neutrality
    
    %initial assumption of activity coefficients
    gCa = 1;   gSi = 1;   gOH = 1;
    
    % thus initial activities and their product
    aCa = gCa*Ca;     aSi = gSi*Si;    aOH = gOH*OH;
    Q = (aCa^3)*(aSi)*(aOH^5);  % having set Si to zero, necessariliy leads to Q = 0
    
    %   iterative loop
    tol =  0.01;  % tolerance for convergence
    count = 0;
    while abs (Kc3s - Q)/Kc3s > tol
        
        count = count + 1;
        % update concentrations of Si and OH to get Q closer to Kc3s
        Qnew = Q + (Kc3s - Q)/100;
        if count == 100000
            Si = -1;
            OH = -1;
            break
        end
        Si = Qnew/(aCa^3)/(aOH^5)/gSi;
        OH = 2*Ca-Si;   % in general, OH could be function of Si, in which case it needs to be updated here   
        
        % compute activity coefficients g via Debye Huckel
        I = 1/2*(Ca*zCa*zCa + Si*zSi*zSi + OH*zOH*zOH); %ionic strength
        gCa = 10^(-A*zCa*zCa*sqrt(I)/(1+B*dCa*sqrt(I))); %coefficients
        gSi = 10^(-A*zSi*zSi*sqrt(I)/(1+B*dSi*sqrt(I)));
        gOH = 10^(-A*zOH*zOH*sqrt(I)/(1+B*dOH*sqrt(I)));
        
        % compute activities and their product
        aCa = gCa * Ca;
        aSi = gSi * Si;
        aOH = gOH * OH;
        Q = (aCa^3)*(aSi)*(aOH^5);
    end
    vSi(i) = Si;
    vOH(i) = OH;
end


semilogy(vCa*1e3,vSi*1e6,'LineWidth',2)
hold on









% CSHI with H3SiO4 ---   CSHI = (Ca + H3SiO4 + OH + 4H2O) with lnK = -22.62     

% set target equilibrium constants
KcshI = exp(-22.62);  %from Bullard 2013

% charges
zCa = 2;
zSi = -1;
zOH = -1;

% hydrated radii of molecules
dCa = 0.82;   %nm
dSi = 0.983;   %computed for H4SiO4, not for H3SiO4, but I guess they are similar enough
dOH = 0.416;


% create vector of molar concentrations in solution
vCa = linspace(0.1,25,100)*1e-3;  % this is fixed; the others computed
vSi = linspace(-1,-1,100);

% for each Ca concentration in vector
for i=1:length(vCa)
    Ca = vCa(i);
    %vtSi = linspace(1e-11,2*vCa(i),1000000);
    vtSi = logspace(-11,log10(2*vCa(i)),1000);
    above = false;
    j=1;
    % scan through the possible values of Si
    while j<length(vtSi)+1 && above == false
        Si = vtSi(j);
        % compute OH as function of Ca and Si
        OH = 2*Ca - Si;
        
         % compute activity coefficients g via Debye Huckel
        I = 1/2*(Ca*zCa*zCa + Si*zSi*zSi + OH*zOH*zOH); %ionic strength
        gCa = 10^(-A*zCa*zCa*sqrt(I)/(1+B*dCa*sqrt(I))); %coefficients
        gSi = 10^(-A*zSi*zSi*sqrt(I)/(1+B*dSi*sqrt(I)));
        gOH = 10^(-A*zOH*zOH*sqrt(I)/(1+B*dOH*sqrt(I)));
        
        % compute activities and their product
        aCa = gCa * Ca;
        aSi = gSi * Si;
        aOH = gOH * OH;
        Q = (aCa)*(aSi)*(aOH);
        
        if Q > KcshI
            above = true;
            vSi(i) = Si;
        end
        j = j + 1;
    end
end

semilogy(vCa*1e3,vSi*1e6)
hold on




%{
vSi = linspace(-1,-1,100);
vOH = linspace(-1,-1,100);

QnewV = zeros(1000,1);

% for each Ca concentration in vector
for i=1:length(vCa)
    
    Ca = vCa(i);    
    Si = 0;  % initial assumption of no Si
    OH = 2*Ca-Si; % OH is linked to Ca and Si content to keep charge neutrality
    
    %initial assumption of activity coefficients
    gCa = 1;   gSi = 1;   gOH = 1;
    
    % thus initial activities and their product
    aCa = gCa*Ca;     aSi = gSi*Si;    aOH = gOH*OH;
    Q = (aCa)*(aSi)*(aOH);  % having set Si to zero, necessariliy leads to Q = 0
    
    %   iterative loop
    tol =  0.01;  % tolerance for convergence
    count = 0;
    while abs (KcshI - Q)/KcshI > tol
        
        count = count + 1;
        % update concentrations of Si and OH to get Q closer to Kc3s
        Qnew = Q + (Kc3s - Q)/2;
        if i== 50
            QnewV(count) = Qnew;
        end
        if count == 1000
            Qnew
            Ca 
            gCa
            Si
            gSi
            gOH
            Si = -1;
            OH = -1;
            break
        end
        Si = Qnew/(aCa)/(aOH)/gSi;
        OH = 2*Ca-Si;   % in general, OH could be function of Si, in which case it needs to be updated here   
        
        % compute activity coefficients g via Debye Huckel
        I = 1/2*(Ca*zCa*zCa + Si*zSi*zSi + OH*zOH*zOH); %ionic strength
        gCa = 10^(-A*zCa*zCa*sqrt(I)/(1+B*dCa*sqrt(I))); %coefficients
        gSi = 10^(-A*zSi*zSi*sqrt(I)/(1+B*dSi*sqrt(I)));
        gOH = 10^(-A*zOH*zOH*sqrt(I)/(1+B*dOH*sqrt(I)));
        
        % compute activities and their product
        aCa = gCa * Ca;
        aSi = gSi * Si;
        aOH = gOH * OH;
        Q = (aCa)*(aSi)*(aOH);
    end
    vSi(i) = Si;
    vOH(i) = OH;
end


semilogy(vCa*1e3,vSi*1e6)
hold on

%}



% CSHII with H3SiO4 --- CSHII = 2Ca + H3SiO4 + 3OH  +3H2O   with lnK = -34.54
% set target equilibrium constants
KcshII = exp(-34.54);  %from Bullard 2013

% charges
zCa = 2;
zSi = -1;
zOH = -1;

% hydrated radii of molecules
dCa = 0.82;   %nm
dSi = 0.983;   %computed for H4SiO4, not for H3SiO4, but I guess they are similar enough
dOH = 0.416;



% create vector of molar concentrations in solution
vCa = linspace(0.1,25,100)*1e-3;  % this is fixed; the others computed
vSi = linspace(-1,-1,100);


% for each Ca concentration in vector
for i=1:length(vCa)
    Ca = vCa(i);
    %vtSi = linspace(1e-11,2*vCa(i),1000000);
    vtSi = logspace(-11,log10(2*vCa(i)),1000);
    above = false;
    j=1;
    % scan through the possible values of Si
    while j<length(vtSi)+1 && above == false
        Si = vtSi(j);
        % compute OH as function of Ca and Si
        OH = 2*Ca - Si;
        
         % compute activity coefficients g via Debye Huckel
        I = 1/2*(Ca*zCa*zCa + Si*zSi*zSi + OH*zOH*zOH); %ionic strength
        gCa = 10^(-A*zCa*zCa*sqrt(I)/(1+B*dCa*sqrt(I))); %coefficients
        gSi = 10^(-A*zSi*zSi*sqrt(I)/(1+B*dSi*sqrt(I)));
        gOH = 10^(-A*zOH*zOH*sqrt(I)/(1+B*dOH*sqrt(I)));
        
        % compute activities and their product
        aCa = gCa * Ca;
        aSi = gSi * Si;
        aOH = gOH * OH;
        Q = (aCa^2)*(aSi)*(aOH^3);
        
        if Q > KcshII
            above = true;
            vSi(i) = Si;
        end
        j = j + 1;
    end
end

semilogy(vCa*1e3,vSi*1e6)
%hold on




%{
% create vector of molar concentrations in solution
vCa = linspace(0.1,25,100)*1e-3;  % this is fixed; the others computed
vSi = linspace(-1,-1,100);
vOH = linspace(-1,-1,100);

% for each Ca concentration in vector
for i=1:length(vCa)
    Ca = vCa(i);
    Si = 0;  % initial assumption of no Si
    OH = 2*Ca-Si; % OH is linked to Ca and Si content to keep charge neutrality
    
    %initial assumption of activity coefficients
    gCa = 1;   gSi = 1;   gOH = 1;
    
    % thus initial activities and their product
    aCa = gCa*Ca;     aSi = gSi*Si;    aOH = gOH*OH;
    Q = (aCa^2)*(aSi)*(aOH^3);  % having set Si to zero, necessariliy leads to Q = 0
    
    %   iterative loop
    tol =  0.01;  % tolerance for convergence
    count = 0;
    while abs (KcshII - Q)/KcshII > tol
        
        count = count + 1 ;
        
        % update concentrations of Si and OH to get Q closer to Kc3s
         Qnew = Q + (Kc3s - Q)/100;
        if count == 100000
            Si = -1;
            OH = -1;
            break
        end
        Si = Qnew/(aCa^2)/(aOH^3)/gSi;
        OH = 2*Ca-Si;   % in general, OH could be function of Si, in which case it needs to be updated here   
        
        % compute activity coefficients g via Debye Huckel
        I = 1/2*(Ca*zCa*zCa + Si*zSi*zSi + OH*zOH*zOH); %ionic strength
        gCa = 10^(-A*zCa*zCa*sqrt(I)/(1+B*dCa*sqrt(I))); %coefficients
        gSi = 10^(-A*zSi*zSi*sqrt(I)/(1+B*dSi*sqrt(I)));
        gOH = 10^(-A*zOH*zOH*sqrt(I)/(1+B*dOH*sqrt(I)));
        
        % compute activities and their product
        aCa = gCa * Ca;
        aSi = gSi * Si;
        aOH = gOH * OH;
        Q = (aCa^2)*(aSi)*(aOH^3);
    end
    vSi(i) = Si;
    vOH(i) = OH;
end


semilogy(vCa*1e3,vSi*1e6)
%}

legend({'Jennings1986','C3S with K = -50.7','C3S with K = -63','CSHI','CSHII'},'FontSize',14)
title('Solubilities assuming only H3SiO4 in solution','FontSize',18)
xlabel('Ca (mmol)', 'FontSize', 18);
ylabel('Si (umol)', 'FontSize', 18);
set(gca,'FontSize',14)
xlim([0 27])
ylim([1e-5 1e10])















% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************

% ***********************************************************************
% ***********************************************************************





% for selected Ca concentrations, plot Keq vs Si concentrations

subplot(1,3,3);



% CSHI with H3SiO4 ---   CSHI = (Ca + H3SiO4 + OH + 4H2O) with lnK = -22.62     

% set target equilibrium constants
KcshI = exp(-22.62);  %from Bullard 2013

% charges
zCa = 2;
zSi = -1;
zOH = -1;

% hydrated radii of molecules
dCa = 0.82;   %nm
dSi = 0.983;   %computed for H4SiO4, not for H3SiO4, but I guess they are similar enough
dOH = 0.416;

% create vector of molar concentrations in solution
vCa = linspace(0.1,25,10)*1e-3;  % one curve will be drawn for each Ca concentration here
vQ = linspace(-1,-1,length(vSi));

% for each Ca concentration in vector
for i=1:length(vCa)
    Ca = vCa(i);
    vSi = linspace(1e-16,2*vCa(i),1000);
    % scan through the possible values of Si
    for j=1:length(vSi)
        Si = vSi(j);
        % compute OH as function of Ca and Si
        OH = 2*Ca - Si;
        
         % compute activity coefficients g via Debye Huckel
        I = 1/2*(Ca*zCa*zCa + Si*zSi*zSi + OH*zOH*zOH); %ionic strength
        gCa = 10^(-A*zCa*zCa*sqrt(I)/(1+B*dCa*sqrt(I))); %coefficients
        gSi = 10^(-A*zSi*zSi*sqrt(I)/(1+B*dSi*sqrt(I)));
        gOH = 10^(-A*zOH*zOH*sqrt(I)/(1+B*dOH*sqrt(I)));
        
        % compute activities and their product
        aCa = gCa * Ca;
        aSi = gSi * Si;
        aOH = gOH * OH;
        vQ(j) = (aCa)*(aSi)*(aOH);
    end
    
    loglog(vSi*1e6,vQ);
    hold on
end
legend({'Ca','C3S with K = -50.7','C3S with K = -63','CSHI','CSHII'},'FontSize',14)


legend({'Ca 1mmol','2.9','5.6','11.2','13.9','16.7','19.5','22.2','25'},'FontSize',14)
title('CSHI activity with H3SiO4 only','FontSize',18)
xlabel('Si (umol)', 'FontSize', 18);
ylabel('Q', 'FontSize', 18);
set(gca,'FontSize',14)
%xlim([0 27])
%ylim([1e-5 1e10])
