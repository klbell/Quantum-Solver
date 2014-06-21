%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mystical Band Diagram Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% By: Kevan Bell
% Last Modified: May 20, 2014
%
%

function BandProgramFinal3()
clc; clear; % clear console and variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nBaseMin            =       0; % Minimum Base Sweep Voltage
nBaseMax            =       1.5; % Maximum Base Sweep Voltage
nCollectorMin       =       1; % Minimum Collector sweep Voltage
nCollectorMax       =       3; % Maximum Collector sweep Voltage
nCollectorStep      =       0.2; % Steps between collector voltages
nStepMin            =       0; % Stepped variable min
nStepMax            =       0; % Stepped varriable max
nStepStep           =       1; % Divitions between steps
WhatStep            =       'T'; %E,C,B = emitter,colector,base thicknesses
sName               =       'Hf/ZrO2/Pt/ZrO2/Al2O3/Hf/300K'; % Gen title

GraphQuality    =       200; % Points per plot (good is around 500)
SimQuality      =       500; % Number of potnetial divisions (around 1000)
yMax            =       4000; % plot y max
T               =       300; % Temp (K)
V1              =       0; % Set emitter voltage
xSpace          =       linspace(nBaseMin,nBaseMax,GraphQuality);

linS = {'-k','--k',':k','-.k','-b','--b',':b','-.b',...
    '-r','--r',':r','-.r','-g','--g',':g','-.g',...
    '-m','--m',':m','-.m'}; % line style list for plots

nCount          =       1;
nColCount       =       1;
linSCount       =       1;
nName           =       1;
bShowPots       =       false; % Display gen potential diagram
bShowU          =       true; % Display U(x)
fignum          =       3; % Figure number counter
figCount        =       1; % Fig counter for filename

for nStep = nStepMin:nStepStep:nStepMax
    for V3 = nCollectorMin:nCollectorStep:nCollectorMax
        for V2 = linspace(nBaseMin,nBaseMax,GraphQuality)
            
            [J(nCount,nColCount),bShowPots,bShowU,LayerThicknesses] =...
                Main(V1,V2,V3,SimQuality,bShowPots,bShowU,nStep,T);
            nCount = nCount+1;
            
        end
        
        Names{nName} = strcat(num2str(V3),' V');
        nName = nName+1;
        
        h = figure(fignum);
        hold on;
        plot(xSpace,J(:,nColCount),linS{linSCount},'LineWidth',1.5);
        
        nCount = 1;
        linSCount = linSCount+1;
        
        if size(linS,2) < linSCount
            
            linSCount = 1; % reset color
            
        end
        nColCount = nColCount+1;
    end
    
    
    switch WhatStep
    
        case 'B' % Base thickness sweep
            
            stitle = strcat(sName,' - 3/15/15 nm, Base: ',...
                num2str(nStep)','nm');
            
            sFileName = strcat('IV 300K',num2str(figCount),'_3_15_15 Base',...
                num2str(nStep),'.png');
            
        case 'C' % Collector thickness sweep
            
            stitle = strcat(sName,' - 3','/',...
                num2str(nStep),'/',num2str(nStep),', Base: 10nm');
            
            sFileName = strcat('IV 300K',num2str(figCount),'_3',...
                '_',num2str(nStep),'_',num2str(nStep));
            
        case 'E' % Emiter thickness sweep
            
            stitle = strcat(sName,' - ',num2str(nStep),'/',...
                '3, Base: 10nm');
            
            sFileName = strcat('IV 300K',num2str(figCount),...
                '_',num2str(nStep),'_15_15');
            
        case 'T' % Temp sweep
            
            stitle = strcat(sName,'-3/3, Base: 10nm ', num2str(T),'K');
            
            sFileName = strcat('IV ',num2str(T),'K - 3_15_15');            
            
    end
    
    sFullFileName = strcat(sFileName,'.png');
    
    legend('Location','NorthWest',Names);
    title(stitle);
    xlabel('Base Voltage (V)');
    ylabel('Current Density (A/cm^2)')
    grid on;
    axis([0 V2 0 yMax]);
    print('-dpng',sFullFileName,'-r300');
%     save(strcat(sFileName,'.mat'),'J');
    close(h);
    
    fignum = fignum+1;
    figCount = figCount+1;
    linSCount = 1; % reset color 
    nColCount = 1; % reset collector counter
end
end


function [J,bShowPot,bShowU,LayerThicknesses] =...
    Main(V1,V2,V3,N,bShowPot,bShowU,nStep,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Setup = struct(...% Input starting information here.
    'L1', struct('Mat',     'Hf',   'AppVolt',  V1   ),...    
    'L2', struct('Mat',     'ZrO2'                  ),...    
    'L3', struct('Mat',     'Pt',   'AppVolt',  V2   ),...
    'L4', struct('Mat',     'ZrO2'                  ),... 
    'L5', struct('Mat',     'Al2O3'                  ),...
    'L6', struct('Mat',     'Hf',   'AppVolt',  V3   ));
    

LayerThicknesses = [2 3 10 1.5 1.5 2].*10^-9; % in nm

Work = struct(... % List of metal work functions (eV).
    'Al', 4.2,...
    'Hf', 3.9,...
    'Co', 5,...
    'Ni', 5.15,...
    'Mn', 4.1,...
    'Mo', 4.65,...
    'Nd', 4.45,...
    'Cr', 4.5,...
    'Ta', 4.4,...
    'Pd', 5.4,...
    'Pt', 5.3,...
    'W',  4.8,...
    'Ti', 4.33,...
    'Ru', 4.71,...
    'V',  4.3,...
    'Zr', 4.05);

Aff = struct(... % List of insulator electron affinities (eV).
    'Al2O3',    1.35,...
    'AlN',      1.9,...
    'HfO2',     2.65,...
    'ZrO2',     2.75,...
    'ZnO',      4.6);

Cap = struct(... % List of insulator capancitances (F/cm^2).
    'Al2O3',    3.98e-6,...
    'AlN',      3.874e-6,...
    'HfO2',     1.107e-5,...
    'ZrO2',     1.107e-5,...
    'ZnO',      3.874e-6);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract Info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[aWorks,aAff,aLogic,nCond,nIns,aVolt,aThick,aSlope] =...
    ExtractInfo(Setup,Work,Aff,Cap);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xMetals,yMetals,xIns,yIns] =...
    PlotBands(nCond,nIns,aLogic,aWorks,aVolt,aAff,aThick,aSlope,Setup,...
    bShowPot);

if bShowPot
    bShowPot = false; % only display static potential
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine Current
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

voltsize = size(aVolt,2);
nVolt = aVolt(voltsize)-aVolt(1);
[J,bShowU] =...
    CurrSolve(LayerThicknesses,nVolt,aVolt,yMetals,yIns,aLogic,N,bShowU,T);

end % end program


function [aWorks,aAff,aLogic,nCond,nIns,aVolt,aThick,aSlope] =...
    ExtractInfo(Setup,Work,Aff,Cap)


nLevel      =       1       ; % Start level counter at 1
nCond       =       1       ; % Start Conductor counter
nIns        =       1       ; % Start Insulator counter
bComplete   =       false   ; % Not finished

nTotLevels  =       size(fieldnames(Setup),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic Information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while ~bComplete
    
    % Find Material, save information
    sName = strcat('Setup.L',num2str(nLevel),'.Mat');
    
    try
        % is metal?
        sTry = eval(strcat('Setup.L',num2str(nLevel),'.AppVolt'));
        bMetal = true;
    catch
        bMetal = false;
    end
    
    if bMetal % if metal
        
        sMetal          =   eval(sName); % metal name
        aWorks(nCond)   =   eval(strcat('Work.',sMetal)); % Save work func
        aLogic(nLevel)  =   1; % Logic 1 means metal layer
        aVolt(nCond)    =   eval(strcat('Setup.L',num2str(nLevel),...
            '.AppVolt')); % applied votages
        
        nCond           =   nCond + 1; % Increment conductor counter
        
    else % if insulator
        
        sIns            =   eval(sName); % insulator name
        aAff(nIns)      =   eval(strcat('Aff.',sIns)); % Save electron aff
        aCap(nIns)      =   eval(strcat('Cap.',sIns)); % Save capacitance
        aLogic(nLevel)  =   0; % Logic 0 means insulator layer
        
        nIns            =   nIns + 1; % Increment insulator counter
        
    end % end if metal
    
    if nLevel == nTotLevels % if at last level
        
        bComplete = true; % now complete
        
    end
    
    nLevel = nLevel + 1; % Increment level
    
end % End while not complete


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thickness of insulator layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nCountThick = 1; % Set counters
nCountLayer = 1;

for i = 1:size(aLogic,2)-1
    
    if aLogic(i) == 0 && aLogic(i+1) == 0 % if current, and next is ins
        
        nCountThick             =   nCountThick + 1; % add to count
        
    elseif aLogic(i) == 0 % if current but not is ins
        
        aThick(nCountLayer)     =   nCountThick; % paste current count
        nCountThick             =   1; % Reset counter
        nCountLayer             =   nCountLayer + 1; % increment ins
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Insulator Slopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aGenSlope = zeros(size(aVolt,1)); % generate empty matrix

for n = 1:size(aVolt,2)-1
    
    % Find slope if capacitances were same, or single level
    aGenSlope(n) = (((-aVolt(n+1)-aWorks(n+1))-(-aVolt(n)-aWorks(n)))/...
        aThick(n))-(aWorks(n)-aWorks(n+1));
    
end

i = 1; % Insulator counter
m = 1; % Metal Counter

for n = 1:size(aThick,2)
    if aThick(n) == 1 % if single layer, then slope is the gen slope
        
        aSlope(i)   =   aGenSlope(n);
        i           =   i + 1; % Increment
        m           =   m + 1; % Increment
        
    elseif aThick(n) == 2 % if double layer
        
        nRatio      =   aCap(i+1)/aCap(i); % Ratio of capacitances
        nMidVoltage =   (aVolt(m+1)-aVolt(m))*(nRatio/(1+nRatio))+aVolt(m);
        % Mean voltage between two insulator layers
        
        if nMidVoltage < aVolt(m) % Classify and apply correct local slope
            
            if nMidVoltage > 0
                
                aSlope(i)   =   aVolt(m)-nMidVoltage;
                aSlope(i+1) =   nMidVoltage-aVolt(m+1);
                
            else
                
                aSlope(i)   =   -(nMidVoltage-aVolt(m));
                aSlope(i+1) =   nMidVoltage-aVolt(m+1);
                
            end
            
        elseif nMidVoltage > aVolt(m)
            
            if nMidVoltage > 0
                
                aSlope(i)   =   -(nMidVoltage-aVolt(m));
                aSlope(i+1) =   nMidVoltage-aVolt(m+1);
                
            else
                
                aSlope(i)   =   aVolt(m)-nMidVoltage;
                aSlope(i+1) =   nMidVoltage-aVolt(m+1);
                
            end
            
        else % if slope is zero, set to zero
            
            aSlope(i)   = 0;
            aSlope(i+1) = 0;
            
        end
        
        i = i + 2; % Increment ins counter
        m = m + 1; % Increment metal counter
        
    end % end of thickness
end % end of layer


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corrections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nCond   =   nCond - 1; % Offset total conductor layers
nIns    =   nIns - 1; % Offset total insulator layers

end % End of gather information

function [xMetals,yMetals,xIns,yIns] =...
    PlotBands(nCond,nIns,aLogic,aWorks,aVolt,aAff,aThick,aSlope,Setup,...
    bShowPot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate X Data For Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xMetals     =   zeros(nCond,2); % Generate empty x metals for plot
xIns        =   zeros(nIns,4); % Generate empty x insulators for plot

m           =   1; % Metal counter
i           =   1; % Insulator counter

for n = 1:size(aLogic,2)
    if aLogic(n) == 1 % If metal
        
        xMetals(m,:)    =   [n-1 n] ; % Set x positions of metals
        m               =   m + 1   ; % Increment metal counter
        
    else % If insulator
        
        xIns(i,:)       =   [n-1 n-1 n n]   ; % Set x of insulators
        i               =   i + 1           ; % Increment ins counter
        
    end % End if metal
end % end for


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Y Data For Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yMetals     =   zeros(nCond,2); % Generate empty y metals for plot
yIns        =   zeros(nIns,4); % Generate empty y insulators for plot

m           =   1; % Metal counter
i           =   1; % Insulator counter
ii          =   1; % Insulator layer counter
bFirst      =   true;

for n = 1:size(aLogic,2)
    if aLogic(n) == 1 % If metal
        
        % Apply driving voltage and work function
        yMetals(m,:)    =   [-aWorks(m)-aVolt(m) -aWorks(m)-aVolt(m)];
        m               =   m + 1; % Increment metal counter
        
    else % if insulator
        if aThick(ii) == 1 % single layer
            
            % Solve y positions
            y1          =   yMetals(m-1,2);
            y2          =   yMetals(m-1,2)+abs(aWorks(m-1)-aAff(i));
            y3          =   y2 + aSlope(i);
            y4          =   y3 - abs(aWorks(m)-aAff(i));
            
            yIns(i,:)   =   [y1 y2 y3 y4]; % Build y array
            
            i           =   i + 1; % Increment ins
            
        else % Multilayer
            if bFirst % if first layer of multilayer
                
                nSaveThick  =   aThick(ii); % Save num of ins layers
                bFirst      =   false; % no longer first layer
                nLayerCount =   1; % Start layer counter
                
                % Solve y positions
                y1          =   yMetals(m-1,2);
                y2          =   yMetals(m-1,2)+abs(aWorks(m-1)-aAff(i));
                y3          =   y2 + aSlope(i);
                
                if aAff(i+1) > aAff(i)
                    y4      =   y3 - abs(aAff(i)-aAff(i+1));
                else
                    y4      =   y3;
                end
                
                nLayerCount     =   nLayerCount + 1; % increment
                
            else % if not the first layer of multilayer
                
                if nLayerCount == nSaveThick % if at last insulator layer
                    
                    bFirst = true; % reset
                    ii = ii + 1; % increment insulator section counter
                    
                end % End if at last layer
                
                % Solve y positions
                y1          =   yIns(i-1,3);
                
                if aAff(i) < aAff(i-1)
                    y2      =   y1+abs(aAff(i)-aAff(i-1));
                else
                    y2      =   y1-abs(aAff(i)-aAff(i-1));
                end
                y3          =   y2 + aSlope(i);
                y4          =   y3 - abs(aAff(i)-aWorks(m));
                
            end % End of if first layer
            
            yIns(i,:) = [y1 y2 y3 y4]; % Build y array
            
            i = i + 1; % Increment
            
        end % end of is single layer?
    end % end of is metal?
end % end of layers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bShowPot
    figure(1)
    hold on; % hold plot on
    for m = 1:size(aVolt,2)
        
        plot(xMetals(m,:),yMetals(m,:),'r','LineWidth',2); % Plot metals
        
    end
    
    for i = 1:size(aAff,2)
        
        plot(xIns(i,:),yIns(i,:),'k','LineWidth',2); % Plot insulators
        
    end
    
    grid on; % Turn on plot grid
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mCount = 1;
iCount = 1;

if bShowPot
    for n = 1:size(fieldnames(Setup))
        
        name = eval(strcat('Setup.L',num2str(n),'.Mat')); % name material
        
        if aLogic(n) == 1 % if metal
            
            text(xMetals(mCount,1),yMetals(mCount,1)-(0.05*min(min(yMetals))),...
                sprintf('%6s',name),'FontSize',12); % Print material name
            
            message = strcat(num2str(eval(strcat('Setup.L',num2str(n),...
                '.AppVolt'))),' V'); % name voltage
            
            text(xMetals(mCount,1),yMetals(mCount,1)-(0.10*min(min(yMetals))),...
                sprintf('%7s',message),'FontSize',12); % print voltage
            
            mCount = mCount + 1; % increment
            
        else
            
            text(xIns(iCount,2),min(yIns(iCount,2:3))+(0.05*min(min(yIns))),...
                sprintf('%7s',name),'FontSize',12); % print material name
            
            iCount = iCount + 1; % increment
            
        end % end if metal
    end % end layers
end
end

function [J,bShowU] =...
    CurrSolve(aThick,nVolt,aVolt,yMetals,yIns,aLogic,N,bShowU,T)

nLvls = size(aThick,2);
nTotLvls = nLvls - 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants and initial values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hbar        =       1.054e-34; % Planks in Js
m           =       9.1094e-31; % rest mass of electron
q           =       1.6e-19; % elementry charge
kB          =       1.38e-23; % Bolztmanns
E           =       0.05*q; % Electron energy
V           =       nVolt*q; % Applied voltage
x           =       linspace(0,sum(aThick(2:nLvls-1)),N+2); % Generate x array
meff        =       m*(0.067 + 0.083*x); % Generate effective mass array

% f=@(x) (3.1*q) - (V/nThick).*x; % square, 0.375eV tall

U = zeros(N+2,1);
k = zeros(N+2,1);
S = zeros(N+1,1);
M = zeros(2,2,N);

% Generate U array

met = 1;
ins = 1;
nCount = 2;
bFirst = true;
metsize = size(yMetals,1);

aThicknm = aThick/10^-9;
NsPerDiv = round(N / sum(aThicknm(2:nTotLvls+1)));

for n = 2:nTotLvls+1
NsPerLvl(n-1) = round(aThicknm(n) * NsPerDiv)+1;
end

CurrLvl = 1;
nLocation = 2;

for n = 1:nTotLvls   
    if aLogic(nCount) == 0 % if insulator
        if bFirst % if first insulator layer
            
            U(1) = yMetals(met,1)*q;
            nSlope = (yIns(ins,3) - yIns(ins,2))*q/aThick(nCount);
            
            if ins ~= size(yIns,1)
                
                for n = 2:NsPerLvl(CurrLvl)+1
                    U(n) = nSlope * x(n) + yIns(ins,2)*q;
                end
                
                nLocation = nLocation + NsPerLvl(CurrLvl);
                
            else
                for n = 2:size(x,2)
                    U(n) = nSlope * x(n) + yIns(ins,2)*q;
                end
            end
            
            U(N+2) = yMetals(metsize,1)*q;
            
            nCount = nCount + 1;
            ins = ins+1;
            bFirst = false;
            CurrLvl = CurrLvl + 1;
            
        else % if not first insulator layer
            
            nSlope = (yIns(ins,3) - yIns(ins,2))*q/aThick(nCount);
            b = true;
            if ins ~= size(yIns,1)
                
                for n = nLocation:nLocation + NsPerLvl(CurrLvl)
                    if b
                        nsave = n;
                        b = false;
                    end
                    U(n) = nSlope * (x(n)-x(nsave)) + yIns(ins,2)*q;
                end
                
            else
                
%                 x(n)-x((nCount-1)*NsPerLvl(CurrLvl))
                for n = nLocation:size(x,2)
                    if b
                        nsave = n;
                        b = false;
                    end
                    U(n) = nSlope * (x(n)-x(nsave)) + yIns(ins,2)*q;
                end
            end
            
            ins = ins+1;
            nLocation = nLocation + NsPerLvl(CurrLvl)+1;
            nCount = nCount+1;
            CurrLvl = CurrLvl + 1;
        end
        
    elseif aLogic(nCount) == 1 % if metal
        met = met+1;
        for n = nLocation:nLocation+NsPerLvl(CurrLvl);
           U(n) = yMetals(met,1)*q;
        end    
        
        nLocation = nLocation + NsPerLvl(CurrLvl); 
        CurrLvl = CurrLvl + 1;      
        nCount = nCount+1; 
    end
   
end

U(n) = yMetals(met+1,1)*q;
U = U-U(1);

% Generate k array
for n = 1:N+2
    k(n) = sqrt(2*meff(n)*(E-U(n)))/hbar;
end

% Generate S array
for n = 1:N+1
    S(n) = (meff(n+1)/meff(n))*(k(n)/k(n+1));
end

% Fill M matries
for n = 1:N+1
    a = (1+S(n))*exp(-i*(k(n+1)-k(n))*x(n));
    b = (1-S(n))*exp(-i*(k(n+1)+k(n))*x(n));
    c = (1-S(n))*exp(i*(k(n+1)+k(n))*x(n));
    d = (1+S(n))*exp(i*(k(n+1)-k(n))*x(n));   
        
    M(:,:,n) = 0.5 * [a b; c d];       
end

% Combine M matrices
Mt = M(:,:,N+1);
for n = 1:N
    Mt = Mt*M(:,:,N+1-n);
end

AN1 = (meff(N+2) * k(1))/(meff(1) * k(N+2) * Mt(2,2));
DE = ((meff(1) * k(N+2))/(meff(N+2) * k(1)))*abs(AN1)^2;
 

J = (((q^2*meff(1)*kB*T)/(2*pi^2*hbar^3))*DE*...
    log((1+exp(-E/(kB*T)))/(1+exp((-V-E)/(kB*T)))))/(100^2);

if bShowU
    figure(2)
    plot(x,U,'k','LineWidth',2);
    title('U(x)');
    xlabel('x');
    ylabel('U(x)');
    bShowU = false;
end

end