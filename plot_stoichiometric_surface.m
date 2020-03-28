%% Plot stoichiometric surface
% Written by: 
% Leonor Guedes da Silva
% LeonorGuedesdaSilva@gmail.com
% 
% Runs in MATLAB R2018a for macOSX
% 
% Package requirements:
% - none
% 
% Input:
%       - Flux distributions
%       - Reconciled data
% 
% Output:
%       
%
%%
clearvars
close all
clc

%% REGRESSION?
regression = 0; % 1 YES, 0 NO

%% ERROR BARS?
errorbars = 0; % 1 YES, 0 NO

%% USE RECONCILED DATA?
RECONCILED_data = 1; % 1 TO USE CORRECTED DATA, 0 TO USE RAW DATA

%% OVERLAP Yagci' model?
Yagci_model = 0; % 1 YES, 0 NO

%% H2Prod or full TCA?
H2_TCA = 0; % 1 YES, 0 NO

%% Import simulation data
fba_results_file = 'S_results_FBA.xls';

NUM = xlsread( fba_results_file , 'AcCoAPrCoAprod_minCO2Prod' ); % 'AcCoAPrCoAprod_maxPHA' or 'normal' or 'withH2' or 'withNOG' or 'fullTCA'
% or 'withPH2MV' or 'H2supplied_maxPHA' or 'GlutamateProd_maxPHA'
CO2_optimum = NUM(:,6); %mol
ATP_optimum = NUM(:,17); %mol
AcCoA_optimum = NUM(:,20); %mol
PrCoA_optimum = NUM(:,21); %mol
Glycolysis_optimum = NUM(:,8); %mol
PDH_optimum = NUM(:,10); %mol
TCA_optimum = NUM(:,12); %mol
oxTCA_optimum = NUM(:,13); %mol
redTCA_optimum = NUM(:,14); %mol
TCAGOX_optimum = NUM(:,15); %mol

Gly = NUM(:,2); Ac = NUM(:,5); %mol
Cfed = 2*Ac + 6*Gly; %Cmol
Gly_Cfed = 6*Gly ./ Cfed;
CO2_Cfed_optimum = CO2_optimum ./ Cfed; %Cmol/Cmol
ATP_Cfed_optimum = ATP_optimum ./ Cfed; %ATPmol/Cmol
ATP_Ac_optimum = ATP_optimum ./ Ac; %ATPmol/mol
AcCoA_Cfed_optimum = 2*AcCoA_optimum ./ Cfed; %Cmol/Cmol
PrCoA_Cfed_optimum = 3*PrCoA_optimum ./ Cfed; %Cmol/Cmol
Glycolysis_Cfed_optimum = 4*Glycolysis_optimum ./ Cfed; %Emol/Cmol
PDH_Cfed_optimum = 2*PDH_optimum ./ Cfed; %Emol/Cmol
TCA_Cfed_optimum = 8*TCA_optimum ./ Cfed; %Emol/Cmol
oxTCA_Cfed_optimum = 4*oxTCA_optimum ./ Cfed; %Emol/Cmol
redTCA_Cfed_optimum = -4*redTCA_optimum ./ Cfed; %Emol/Cmol
TCAGOX_Cfed_optimum = 2*TCAGOX_optimum ./ Cfed; %Emol/Cmol

NUM = xlsread( fba_results_file , 'AcCoAPrCoAprod_maxCO2Prod' ); % 'AcCoAPrCoAprod_minPHA' or'ex3' or 'normal_minPHA' or 'H2supplied_maxPHV'
% or 'GlutamateProd_minPHA'
CO2_minimum = NUM(:,6); %mol
ATP_minimum = NUM(:,17); %mol
AcCoA_minimum = NUM(:,20); %mol
PrCoA_minimum = NUM(:,21); %mol
Glycolysis_minimum = NUM(:,8); %mol
PDH_minimum = NUM(:,10); %mol
TCA_minimum = NUM(:,12); %mol
oxTCA_minimum = NUM(:,13); %mol
redTCA_minimum = NUM(:,14); %mol
TCAGOX_minimum = NUM(:,15); %mol

CO2_Cfed_minimum = CO2_minimum ./ Cfed; %Cmol/Cmol
AcCoA_Cfed_minimum = 2*AcCoA_minimum ./ Cfed; %Cmol/Cmol
ATP_Cfed_minimum = ATP_minimum ./ Cfed; %ATPmol/Cmol
ATP_Ac_minimum = ATP_minimum ./ Ac; %ATPmol/mol
PrCoA_Cfed_minimum = 3*PrCoA_minimum ./ Cfed; %Cmol/Cmol
Glycolysis_Cfed_minimum = 4*Glycolysis_minimum ./ Cfed; %Emol/Cmol
PDH_Cfed_minimum = 2*PDH_minimum ./ Cfed; %Emol/Cmol
TCA_Cfed_minimum = 8*TCA_minimum ./ Cfed; %Emol/Cmol
oxTCA_Cfed_minimum = 4*oxTCA_minimum ./ Cfed; %Emol/Cmol
redTCA_Cfed_minimum = -4*redTCA_minimum ./ Cfed; %Emol/Cmol
TCAGOX_Cfed_minimum = 2*TCAGOX_minimum ./ Cfed; %Emol/Cmol

Gly_Cfed = round(Gly_Cfed, 4) - round(6/18, 4); %place origin in stoichiometric molar relation 1Glyc:6Ac

%% Import extra simulation data
if H2_TCA
    
    % H2 production
    NUM = xlsread( fba_results_file , 'H2Prod_minCO2Prod' );
    CO2_optimum_H2 = NUM(:,6); %mol
    AcCoA_optimum_H2 = NUM(:,20); %mol
    PrCoA_optimum_H2 = NUM(:,21); %mol
    
    Gly_H2 = NUM(:,2); Ac_H2 = NUM(:,5); %mol
    Cfed_H2 = 2*Ac_H2 + 6*Gly_H2; %Cmol
    Gly_Cfed_H2 = 6*Gly_H2 ./ Cfed_H2;
    
    CO2_Cfed_optimum_H2 = CO2_optimum_H2 ./ Cfed_H2; %Cmol/Cmol
    AcCoA_Cfed_optimum_H2 = 2*AcCoA_optimum_H2 ./ Cfed_H2; %Cmol/Cmol
    PrCoA_Cfed_optimum_H2 = 3*PrCoA_optimum_H2 ./ Cfed_H2; %Cmol/Cmol
    
    NUM = xlsread( fba_results_file , 'H2Prod_maxCO2Prod' );
    CO2_minimum_H2 = NUM(:,6); %mol
    AcCoA_minimum_H2 = NUM(:,20); %mol
    PrCoA_minimum_H2 = NUM(:,21); %mol
    
    CO2_Cfed_minimum_H2 = CO2_minimum_H2 ./ Cfed_H2; %Cmol/Cmol
    AcCoA_Cfed_minimum_H2 = 2*AcCoA_minimum_H2 ./ Cfed_H2; %Cmol/Cmol
    PrCoA_Cfed_minimum_H2 = 3*PrCoA_minimum_H2 ./ Cfed_H2; %Cmol/Cmol
    
    % full TCA operation
    NUM = xlsread( fba_results_file , 'fullTCA_minCO2Prod' ); % 'fullTCA_maxPHA'
    CO2_optimum_TCA = NUM(:,6); %mol
    AcCoA_optimum_TCA = NUM(:,20); %mol
    PrCoA_optimum_TCA = NUM(:,21); %mol
    
    Gly_TCA = NUM(:,2); Ac_TCA = NUM(:,5); %mol
    Cfed_TCA = 2*Ac_TCA + 6*Gly_TCA; %Cmol
    Gly_Cfed_TCA = 6*Gly_TCA ./ Cfed_TCA;
    
    CO2_Cfed_optimum_TCA = CO2_optimum_TCA ./ Cfed_TCA; %Cmol/Cmol
    AcCoA_Cfed_optimum_TCA = 2*AcCoA_optimum_TCA ./ Cfed_TCA; %Cmol/Cmol
    PrCoA_Cfed_optimum_TCA = 3*PrCoA_optimum_TCA ./ Cfed_TCA; %Cmol/Cmol
    
    NUM = xlsread( fba_results_file , 'fullTCA_maxCO2Prod' ); % 'fullTCA_minPHA'
    CO2_minimum_TCA = NUM(:,6); %mol
    AcCoA_minimum_TCA = NUM(:,20); %mol
    PrCoA_minimum_TCA = NUM(:,21); %mol
    
    CO2_Cfed_minimum_TCA = CO2_minimum_TCA ./ Cfed_TCA; %Cmol/Cmol
    AcCoA_Cfed_minimum_TCA = 2*AcCoA_minimum_TCA ./ Cfed_TCA; %Cmol/Cmol
    PrCoA_Cfed_minimum_TCA = 3*PrCoA_minimum_TCA ./ Cfed_TCA; %Cmol/Cmol
    
end

%% Import experimental data

literature_data_file = 'tableS5_LiteratureData_reconciled.xlsx';
[NUM,TXT] = xlsread( literature_data_file , 'All' );

Ac_exp = NUM(:,3);
Gly_Ac_exp = NUM(:,4);
Cfed_exp = (Ac_exp + Gly_Ac_exp);
Gly_Cfed_exp = Gly_Ac_exp ./ Cfed_exp;
PHV_Cfed_exp = NUM(:,5) ./ Cfed_exp;
PHB_Cfed_exp = NUM(:,6) ./ Cfed_exp;
PH2MV_Cfed_exp = NUM(:,7) ./ Cfed_exp;
AcCoA_Cfed_exp = PHB_Cfed_exp + 2/5 * PHV_Cfed_exp;
PrCoA_Cfed_exp = PH2MV_Cfed_exp + 3/5 * PHV_Cfed_exp;
CO2_Cfed_exp = -NUM(:,8) ./ Cfed_exp; %calculated based on closing the C balance

if RECONCILED_data
    Ac_exp_corrected =      NUM(:,15);
    Gly_exp_corrected =     NUM(:,16);
    Cfed_exp_corrected =    (Ac_exp_corrected + Gly_exp_corrected);
    PHV_exp_corrected =     NUM(:,17);
    PHB_exp_corrected =     NUM(:,18);
    PH2MV_exp_corrected =   NUM(:,19);
    CO2_exp_corrected =     NUM(:,20);
    
    Gly_Cfed_exp_corrected =    Gly_exp_corrected   ./ Cfed_exp_corrected;
    PHV_Cfed_exp_corrected =    PHV_exp_corrected   ./ Cfed_exp_corrected;
    PHB_Cfed_exp_corrected =    PHB_exp_corrected   ./ Cfed_exp_corrected;
    PH2MV_Cfed_exp_corrected =  PH2MV_exp_corrected ./ Cfed_exp_corrected;
    AcCoA_Cfed_exp_corrected = PHB_Cfed_exp_corrected + 2/5 * PHV_Cfed_exp_corrected;
    PrCoA_Cfed_exp_corrected = PH2MV_Cfed_exp_corrected + 3/5 * PHV_Cfed_exp_corrected;
    CO2_Cfed_exp_corrected =    CO2_exp_corrected   ./ Cfed_exp_corrected;
    
    Cfed_exp =          Cfed_exp_corrected;
    Gly_Cfed_exp =      Gly_Cfed_exp_corrected;
    CO2_Cfed_exp =      CO2_Cfed_exp_corrected;
    PHV_Cfed_exp =      PHV_Cfed_exp_corrected;
    PHB_Cfed_exp =      PHB_Cfed_exp_corrected;
    PH2MV_Cfed_exp =    PH2MV_Cfed_exp_corrected;
    AcCoA_Cfed_exp =    AcCoA_Cfed_exp_corrected;
    PrCoA_Cfed_exp =    PrCoA_Cfed_exp_corrected;
    
    % Error propagation
    Ac_exp_STD =    NUM(:,21);
    Gly_exp_STD =   NUM(:,22);
    PHV_exp_STD =   NUM(:,23);
    PHB_exp_STD =   NUM(:,24);
    PH2MV_exp_STD = NUM(:,25);
    CO2_exp_STD =   NUM(:,26);
    
    Cfed_exp_STD =          (Ac_exp_STD.^2 + Gly_exp_STD.^2).^0.5 ;
    Gly_Cfed_exp_STD =      abs(Gly_Cfed_exp) .* ( (Gly_exp_STD./Gly_exp_corrected).^2 + (Cfed_exp_STD./Cfed_exp).^2 ).^0.5;
    CO2_Cfed_exp_STD =      abs(CO2_Cfed_exp) .* ( (CO2_exp_STD./CO2_exp_corrected).^2 + (Cfed_exp_STD./Cfed_exp).^2 ).^0.5;
    PHV_Cfed_exp_STD =      abs(PHV_Cfed_exp) .* ( (PHV_exp_STD./PHV_exp_corrected).^2 + (Cfed_exp_STD./Cfed_exp).^2 ).^0.5;
    PHB_Cfed_exp_STD =      abs(PHB_Cfed_exp) .* ( (PHB_exp_STD./PHB_exp_corrected).^2 + (Cfed_exp_STD./Cfed_exp).^2 ).^0.5;
    PH2MV_Cfed_exp_STD =    abs(PH2MV_Cfed_exp) .* ( (PH2MV_exp_STD./PH2MV_exp_corrected).^2 + (Cfed_exp_STD./Cfed_exp).^2 ).^0.5;
    AcCoA_Cfed_exp_STD =    ( PHB_exp_STD.^2 + (abs(2/5) .* PHV_exp_STD).^2 ).^0.5;
    PrCoA_Cfed_exp_STD =    ( PH2MV_exp_STD.^2 + (abs(3/5) .* PHV_exp_STD).^2 ).^0.5;
    
end

Gly_Cfed_exp = Gly_Cfed_exp -6/18; %place origin in stoichiometric molar relation 1Glyc:6Ac

Author = TXT(2:end,1);
Year = TXT(2:end,2);

%% Yagci model
if Yagci_model
    %PAO
    f_GLX = [linspace(0,10,20) 100 1000];
    Ac_Yagci_PAO = 2 * (6 + 3 * f_GLX); %Cmol
    Gly_Yagci_PAO = 6 * 1; %Cmol
    Cfed_Yagci_PAO = (Ac_Yagci_PAO + Gly_Yagci_PAO);
    PHB_Yagci_PAO = 4 * 4; %Cmol
    PHV_Yagci_PAO = 5 * f_GLX; %Cmol
    CO2_Yagci_PAO = 1 * (2 + f_GLX); %Cmol
    
    Ac_Cfed_Yagci_PAO = Ac_Yagci_PAO ./ Cfed_Yagci_PAO; %Cmol/Cmol
    Gly_Cfed_Yagci_PAO = Gly_Yagci_PAO ./ Cfed_Yagci_PAO; %Cmol/Cmol
    PHB_Cfed_Yagci_PAO = PHB_Yagci_PAO ./ Cfed_Yagci_PAO; %Cmol/Cmol
    PHV_Cfed_Yagci_PAO = PHV_Yagci_PAO ./ Cfed_Yagci_PAO; %Cmol/Cmol
    CO2_Cfed_Yagci_PAO = CO2_Yagci_PAO ./ Cfed_Yagci_PAO; %Cmol/Cmol
    PH2MV_Cfed_Yagci_PAO = 0; %NOT MODELLED
    AcCoA_Cfed_Yagci_PAO = PHB_Cfed_Yagci_PAO + 2/5 * PHV_Cfed_Yagci_PAO;
    PrCoA_Cfed_Yagci_PAO = PH2MV_Cfed_Yagci_PAO + 3/5 * PHV_Cfed_Yagci_PAO;
    
    Gly_Cfed_Yagci_PAO = Gly_Cfed_Yagci_PAO -6/18; %place origin in stoichiometric molar relation 1Glyc:6Ac
    
    %GAO
    a_GAO = [linspace(0,10,20) 100 1000];
    Ac_Yagci_GAO = 2 * 1; %Cmol
    Gly_Yagci_GAO = 6 * (1/3 + 2/3 * a_GAO); %Cmol
    Cfed_Yagci_GAO = (Ac_Yagci_GAO + Gly_Yagci_GAO);
    PHB_Yagci_GAO = 4 * 2/3; %Cmol
    PHV_Yagci_GAO = 5 * (1/6 + 2/3 * a_GAO); %Cmol
    CO2_Yagci_GAO = 1 * (1/2 + 2/3 * a_GAO); %Cmol
    
    Ac_Cfed_Yagci_GAO = Ac_Yagci_GAO ./ Cfed_Yagci_GAO; %Cmol/Cmol
    Gly_Cfed_Yagci_GAO = Gly_Yagci_GAO ./ Cfed_Yagci_GAO; %Cmol/Cmol
    PHB_Cfed_Yagci_GAO = PHB_Yagci_GAO ./ Cfed_Yagci_GAO; %Cmol/Cmol
    PHV_Cfed_Yagci_GAO = PHV_Yagci_GAO ./ Cfed_Yagci_GAO; %Cmol/Cmol
    CO2_Cfed_Yagci_GAO = CO2_Yagci_GAO ./ Cfed_Yagci_GAO; %Cmol/Cmol
    PH2MV_Cfed_Yagci_GAO = 0; %NOT MODELLED
    AcCoA_Cfed_Yagci_GAO = PHB_Cfed_Yagci_GAO + 2/5 * PHV_Cfed_Yagci_GAO;
    PrCoA_Cfed_Yagci_GAO = PH2MV_Cfed_Yagci_GAO + 3/5 * PHV_Cfed_Yagci_GAO;
    
    Gly_Cfed_Yagci_GAO = Gly_Cfed_Yagci_GAO -6/18; %place origin in stoichiometric molar relation 1Glyc:6Ac
end

%% PLOT CO2 YIELD
figure('Name','CO2','NumberTitle','off'), clf, hold on

theoreticalY_CO2_minimumPHB_Cfed = 1-16/18; %cmolCO2/cmol MINIMUM FOR THIS PHA BASED METABOLISM
plot(Gly_Cfed,theoreticalY_CO2_minimumPHB_Cfed*ones(size(Gly_Cfed)),':k', 'linewidth', 2,'Color',[0.75 0.75 0.75]);

plot(Gly_Cfed(Gly_Cfed<=0),CO2_Cfed_optimum(Gly_Cfed<=0),'--r','Color',[226/256 26/256 26/256], 'linewidth', 2)
plot(Gly_Cfed(Gly_Cfed>=0),CO2_Cfed_optimum(Gly_Cfed>=0),'--g','Color',[165/256 202/256 26/256], 'linewidth', 2)
plot(Gly_Cfed,CO2_Cfed_minimum, ':k', 'linewidth', 2);
cmap = parula;
fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
    [CO2_Cfed_minimum;CO2_Cfed_optimum(end:-1:1);CO2_Cfed_minimum(1)],...
    'y','FaceColor',cmap(end-5,:),'FaceAlpha',0.1,'LineStyle','none');

plot(0,theoreticalY_CO2_minimumPHB_Cfed,'p','Color', [0/256 166/256 214/256],'MarkerSize',14, 'linewidth', 2)

if Yagci_model
    plot(Gly_Cfed_Yagci_PAO,CO2_Cfed_Yagci_PAO,':r','Color',[109/256 23/256 127/256], 'linewidth', 2)
    plot(Gly_Cfed_Yagci_GAO,CO2_Cfed_Yagci_GAO,':g','Color',[0/256 136/256 145/256], 'linewidth', 2)
    
    text(-0.3, 0.152, {'Yagci (2003) PAO model'}, 'fontsize', 14,'rotation', -18,...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom','Color',[109/256 23/256 127/256])
    text(0.39, 0.135, {'Yagci (2003) GAO model'}, 'fontsize', 14,'rotation', +10,...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom','Color', [0/256 136/256 145/256])
end

if H2_TCA
    H2_color = [0/256 136/256 145/256];
    plot(Gly_Cfed,CO2_Cfed_optimum_H2,'--r','Color',H2_color, 'linewidth', 2)
    plot(Gly_Cfed,CO2_Cfed_minimum_H2,':g','Color',H2_color, 'linewidth', 2)
    fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
        [CO2_Cfed_minimum_H2;CO2_Cfed_optimum_H2(end:-1:1);CO2_Cfed_minimum_H2(1)],...
        'y','FaceColor',H2_color,'FaceAlpha',0.1,'LineStyle','none');
    text(-0.3, 0.22, {'with H_2 prod.'}, 'fontsize', 14,'rotation', -0,...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom','Color',H2_color)
    
    TCA_color = [109/256 23/256 127/256];
    plot(Gly_Cfed,CO2_Cfed_optimum_TCA,'--r','Color',TCA_color, 'linewidth', 2)
    plot(Gly_Cfed,CO2_Cfed_minimum_TCA,':g','Color',TCA_color, 'linewidth', 2)
    fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
        [CO2_Cfed_minimum_TCA;CO2_Cfed_optimum_TCA(end:-1:1);CO2_Cfed_minimum_TCA(1)],...
        'y','FaceColor',TCA_color,'FaceAlpha',0.1,'LineStyle','none');
    
    text(-0.3, 0.13, {'with full TCA'}, 'fontsize', 14,'rotation', -0,...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom','Color',TCA_color)
end

% get rid of standard axes decorations
set(gca, 'box', 'off')

% Fix the axes sizes
y_axis_size = 0.26;
x_axis_size = 0.8;
axis([-0.5+0.11 x_axis_size-0.41 0.05 y_axis_size])
yl=ylim; xl=xlim;
ax = gca;
ax.YAxisLocation = 'origin';
ax.ColorOrderIndex = 1;

% the static texts
text(0,    y_axis_size, {'CO_2 produced per C consumed','(Cmol/Cmol)'}, 'rotation', 0, 'fontsize', 14,...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom')

text(xl(2) , yl(1)-0.012, 'Stoichiometric proportion C sources consumed (Cmol)', 'fontsize', 14,...
    'HorizontalAlignment','right', 'VerticalAlignment','top')

text(xl(2), 0.145, {'C) Glycolysis +','redTCA branch'}, 'fontsize', 14,'rotation', +10,...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','Color', [165/256 202/256 26/256])

text(xl(1)+0.1, 0.16, {'A) Glyoxylate Shunt +','Glycolysis'}, 'fontsize', 14,'rotation', -18,...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom','Color',[226/256 26/256 26/256])

text(0.005, theoreticalY_CO2_minimumPHB_Cfed-0.005, {'B) Glycolysis only'}, 'fontsize', 14,'rotation', +0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [0/256 166/256 214/256])

% PLOT experimental data
if errorbars
    errorbar(Gly_Cfed_exp,CO2_Cfed_exp,CO2_Cfed_exp_STD,CO2_Cfed_exp_STD,Gly_Cfed_exp_STD,Gly_Cfed_exp_STD,...
        'Color',[0.8 0.8 0.8],'LineStyle','none','LineWidth',0.2)
end

leg = {strcat(num2str(1), '.  ', Author{1},' (',Year{1}, ')')};
counter = 1;
pos=1;
for ii = 2:length(Author)
    if isequal(Author{ii},Author{ii-1}) && isequal(Year{ii},Year{ii-1})
        leg{pos} = strcat(num2str(counter), '-',num2str(ii), '.  ', Author{ii},' (',Year{ii}, ')');
    else
        pos = pos+1;
        counter = ii;
        leg{pos} = strcat(num2str(ii), '.  ', Author{ii},' (',Year{ii}, ')');
    end
end


enrichment = TXT(2:end,3);
PAO_idx = contains(enrichment,'PAO');
PAO1_idx = contains(enrichment,'PAO I ');
PAO2_idx = contains(enrichment,'PAO II');
GAO_idx = contains(enrichment,'GAO');
EBPR_idx = contains(enrichment,'EBPR');
%unidentified = contains(enrichment,'????');


p1 = plot(Gly_Cfed_exp,CO2_Cfed_exp,'s','Color',[0.5 0.5 0.5],'MarkerSize',6);
p2 = plot(Gly_Cfed_exp(PAO_idx),CO2_Cfed_exp(PAO_idx),'o','MarkerSize',6);
p3 = plot(Gly_Cfed_exp(PAO1_idx),CO2_Cfed_exp(PAO1_idx),'.','MarkerSize',12);
p4 = plot(Gly_Cfed_exp(PAO2_idx),CO2_Cfed_exp(PAO2_idx),'.','MarkerSize',6);
p5 = plot(Gly_Cfed_exp(GAO_idx),CO2_Cfed_exp(GAO_idx),'xk','linewidth',1,'MarkerSize',12);
p6 = plot(Gly_Cfed_exp(EBPR_idx),CO2_Cfed_exp(EBPR_idx),'.','MarkerSize',12);

text(Gly_Cfed_exp,CO2_Cfed_exp,num2str([1:length(CO2_Cfed_exp)]'),...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom')

text( xl(1)+0.04 , y_axis_size , leg , ...
    'HorizontalAlignment','right', 'VerticalAlignment','top',...
    'FontSize',12)

% pretty plot
lgd = legend([p1 p2 p3 p4 p5 p6],'Unidentified','PAO','PAO I', 'PAO II', 'GAO','EBPR sludge',...
    'Location','SouthEast'); legend('boxoff');
lgd.FontSize = 14;
title(lgd,'Enrichment type')

set(gca,'color','white')

yl_main=ylim; xl_main=xlim;

xticks([(0-6/18)	(1/6-6/18)     0	(1/2-6/18)     (2/3-6/18)    (5/6-6/18)])
xticklabels({'0 Gly:1 Ac','^{1}/_{6}Gly:^{5}/_{6}Ac','^{2}/_{6}Gly:^{4}/_{6}Ac','^{3}/_{6}Gly:^{3}/_{6}Ac','^{4}/_{6}Gly:^{2}/_{6}Ac','^{5}/_{6}Gly:^{1}/_{6}Ac'})
ax.FontSize = 12;

% make insert of whole triangle
% create smaller axes in top right, and plot on it
axes('Position',[.705 .75 .2 .2])

p1 = plot(Gly_Cfed,theoreticalY_CO2_minimumPHB_Cfed*ones(size(Gly_Cfed)),':k', 'linewidth', 2,'Color',[0.75 0.75 0.75]);
hold on
p2 = plot(Gly_Cfed,CO2_Cfed_optimum,'--k', 'linewidth', 1);
p3 = plot(Gly_Cfed,CO2_Cfed_minimum, ':k', 'linewidth', 2);
p4 = fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
    [CO2_Cfed_minimum;CO2_Cfed_optimum(end:-1:1);CO2_Cfed_minimum(1)],...
    'y','FaceColor',cmap(end-5,:),'FaceAlpha',0.2,'LineStyle','none');

% get rid of standard axes decorations
set(gca, 'box', 'off')

% Fix the axes sizes
y_axis_size = 1.04;
% x_axis_size = 0.8;
axis([-0.5+0.01 x_axis_size-0.01 0 y_axis_size])
yl=ylim; xl=xlim;
ax = gca;
ax.YAxisLocation = 'origin';
ax.ColorOrderIndex = 1;

hold on
% PLOT experimental data
p = plot(Gly_Cfed_exp,CO2_Cfed_exp,'s','Color',[0.5 0.5 0.5],'MarkerSize',3);
p = plot(Gly_Cfed_exp(PAO_idx),CO2_Cfed_exp(PAO_idx),'o','MarkerSize',3);
p = plot(Gly_Cfed_exp(PAO1_idx),CO2_Cfed_exp(PAO1_idx),'.','MarkerSize',6);
p = plot(Gly_Cfed_exp(PAO2_idx),CO2_Cfed_exp(PAO2_idx),'.','MarkerSize',3);
p = plot(Gly_Cfed_exp(GAO_idx),CO2_Cfed_exp(GAO_idx),'xk','linewidth',0.5,'MarkerSize',6);
p = plot(Gly_Cfed_exp(EBPR_idx),CO2_Cfed_exp(EBPR_idx),'.','MarkerSize',6);

set(gca,'color','white')

p = fill([xl_main(1) xl_main(1) xl_main(2) xl_main(2)],...
    [yl_main yl_main(end:-1:1)],...
    'k','FaceAlpha',0.2,'LineStyle','none');

lgd = legend([p1 p2 p3 p4],'Theoretical min','Optimum','Minimum','Solution space',...
    'Location','NorthEast'); legend('boxoff');
lgd.FontSize = 10;
lgd.Position = [0.865-0.05 0.795+0.05 0.0764 0.0427];

xticks([(0-6/18)	0     (2/3-6/18)    (1-6/18)])
xticklabels({'0 Gly:1 Ac','^{2}/_{6}Gly:^{4}/_{6}Ac','^{4}/_{6}Gly:^{2}/_{6}Ac','1 Gly:0 Ac'})

ax.TickLength = [0.02 0.035];

%% PLOT AcCoA YIELD
figure('Name','AcCoA','NumberTitle','off'), clf, hold on
%set(gcf,'units','normalized','outerposition',[0 0 1 1])

theoreticalY_PHB_Cfed = 16/18; %cmolPHB/cmolAc
plot(Gly_Cfed,theoreticalY_PHB_Cfed*ones(size(Gly_Cfed)),':k', 'linewidth', 2,'Color',[0.75 0.75 0.75])

plot(Gly_Cfed(Gly_Cfed<=0),AcCoA_Cfed_optimum(Gly_Cfed<=0),'--r','Color',[226/256 26/256 26/256], 'linewidth', 2)
plot(Gly_Cfed(Gly_Cfed>=0),AcCoA_Cfed_optimum(Gly_Cfed>=0),'--g','Color',[165/256 202/256 26/256], 'linewidth', 2)
plot(Gly_Cfed,AcCoA_Cfed_minimum, ':k', 'linewidth', 2)
cmap = parula;
fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
    [AcCoA_Cfed_minimum;AcCoA_Cfed_optimum(end:-1:1);AcCoA_Cfed_minimum(1)],...
    'y','FaceColor',cmap(end-5,:),'FaceAlpha',0.1,'LineStyle','none')

plot(0,theoreticalY_PHB_Cfed,'p','Color', [0/256 166/256 214/256],'MarkerSize',14, 'linewidth', 2)

if Yagci_model
    plot(Gly_Cfed_Yagci_PAO,AcCoA_Cfed_Yagci_PAO,':r','Color',[109/256 23/256 127/256], 'linewidth', 2)
    plot(Gly_Cfed_Yagci_GAO,AcCoA_Cfed_Yagci_GAO,':g','Color',[0/256 136/256 145/256], 'linewidth', 2)
    
    text(-0.27, 0.44, {'Yagci (2003) PAO model'}, 'fontsize', 14,'rotation', +49,...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom','Color',[109/256 23/256 127/256])
    text(0.39, 0.565, {'Yagci (2003) GAO model'}, 'fontsize', 14,'rotation', -29,...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom','Color', [0/256 136/256 145/256])
end

if H2_TCA
    H2_color = [0/256 136/256 145/256];
    plot(Gly_Cfed,AcCoA_Cfed_optimum_H2,'--r','Color',H2_color, 'linewidth', 2)
    plot(Gly_Cfed,AcCoA_Cfed_minimum_H2,':g','Color',H2_color, 'linewidth', 2)
    fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
        [AcCoA_Cfed_minimum_H2;AcCoA_Cfed_optimum_H2(end:-1:1);AcCoA_Cfed_minimum_H2(1)],...
        'y','FaceColor',H2_color,'FaceAlpha',0.1,'LineStyle','none');
    text(0.27, 0.7, {'with','H_2', 'prod.'}, 'fontsize', 14,'rotation', -0,...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom','Color',H2_color)
    
    TCA_color = [109/256 23/256 127/256];
    plot(Gly_Cfed,AcCoA_Cfed_optimum_TCA,'--r','Color',TCA_color, 'linewidth', 2)
    plot(Gly_Cfed,AcCoA_Cfed_minimum_TCA,':g','Color',TCA_color, 'linewidth', 2)
    fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
        [AcCoA_Cfed_minimum_TCA;AcCoA_Cfed_optimum_TCA(end:-1:1);AcCoA_Cfed_minimum_TCA(1)],...
        'y','FaceColor',TCA_color,'FaceAlpha',0.1,'LineStyle','none');
    
    text(-0.3, 0.7, {'with full TCA'}, 'fontsize', 14,'rotation', -0,...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom','Color',TCA_color)
end

% get rid of standard axes decorations
set(gca, 'box', 'off')

% Fix the axes sizes
y_axis_size = 1.04;
% x_axis_size = 0.5;
axis([-0.5+0.11 x_axis_size-0.41 0.4 y_axis_size])
yl=ylim; xl=xlim;
ax = gca;
ax.YAxisLocation = 'origin';
ax.ColorOrderIndex = 1;

% the static texts
text(0,    y_axis_size, {'AcCoA* accumulated per C consumed','(Cmol/Cmol)'}, 'rotation', 0, 'fontsize', 14,...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom')

text(xl(2) , yl(1)-0.03, 'Stoichiometric proportion C sources consumed (Cmol)', 'fontsize', 14,...
    'HorizontalAlignment','right', 'VerticalAlignment','top')

text(0.36, yl(1)+0.13, {'C) Glycolysis +','redTCA branch'}, 'fontsize', 14, 'rotation', -29, ...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','Color', [165/256 202/256 26/256])

text(-0.24, yl(1)+0.01, {'A) Glyoxylate Shunt +','Glycolysis'}, 'fontsize', 14, 'rotation', +49, ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom','Color',[226/256 26/256 26/256])

text(0.01, 0.9, {'B) Glycolysis only'}, 'fontsize', 14,'rotation', +0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [0/256 166/256 214/256])

% PLOT experimental data
if errorbars
    errorbar(Gly_Cfed_exp,AcCoA_Cfed_exp,AcCoA_Cfed_exp_STD,AcCoA_Cfed_exp_STD,Gly_Cfed_exp_STD,Gly_Cfed_exp_STD,...
        'Color',[0.8 0.8 0.8],'LineStyle','none','LineWidth',0.2)
end

leg = {strcat(num2str(1), '.  ', Author{1},' (',Year{1}, ')')};
counter = 1;
pos=1;
for ii = 2:length(Author)
    if isequal(Author{ii},Author{ii-1}) && isequal(Year{ii},Year{ii-1})
        leg{pos} = strcat(num2str(counter), '-',num2str(ii), '.  ', Author{ii},' (',Year{ii}, ')');
    else
        pos = pos+1;
        counter = ii;
        leg{pos} = strcat(num2str(ii), '.  ', Author{ii},' (',Year{ii}, ')');
    end
end


enrichment = TXT(2:end,3);
PAO_idx = contains(enrichment,'PAO');
PAO1_idx = contains(enrichment,'PAO I ');
PAO2_idx = contains(enrichment,'PAO II');
GAO_idx = contains(enrichment,'GAO');
EBPR_idx = contains(enrichment,'EBPR');
%unidentified = contains(enrichment,'????');


p1 = plot(Gly_Cfed_exp,AcCoA_Cfed_exp,'s','Color',[0.5 0.5 0.5],'MarkerSize',6);
p2 = plot(Gly_Cfed_exp(PAO_idx),AcCoA_Cfed_exp(PAO_idx),'o','MarkerSize',6);
p3 = plot(Gly_Cfed_exp(PAO1_idx),AcCoA_Cfed_exp(PAO1_idx),'.','MarkerSize',12);
p4 = plot(Gly_Cfed_exp(PAO2_idx),AcCoA_Cfed_exp(PAO2_idx),'.','MarkerSize',6);
p5 = plot(Gly_Cfed_exp(GAO_idx),AcCoA_Cfed_exp(GAO_idx),'xk','linewidth',1,'MarkerSize',12);
p6 = plot(Gly_Cfed_exp(EBPR_idx),AcCoA_Cfed_exp(EBPR_idx),'.','MarkerSize',12);


text(Gly_Cfed_exp,AcCoA_Cfed_exp,num2str([1:length(AcCoA_Cfed_exp)]'),...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom')

text( xl(1)+0.04 , yl(2) , leg , ...
    'HorizontalAlignment','right', 'VerticalAlignment','top',...
    'FontSize',12)

if regression
    % Linear regression (for glycogen LIMITATION)
    x = Gly_Cfed_exp(Gly_Cfed_exp<0); y_meas = AcCoA_Cfed_exp(Gly_Cfed_exp<0);
    A = [x ones(size(x))];
    s_y = AcCoA_Cfed_exp_STD(Gly_Cfed_exp<0);
    S_y = diag( s_y.^2 ); % covariance matrix
    beta = inv( A'*inv(S_y)*A) * A' *inv(S_y) * y_meas;
    
    x_fit = -0.4:0.1:0;
    y_fit = [x_fit' ones(size(x_fit'))] * beta;
    plot(x_fit,y_fit,'Color',autumn(1));
    
    r = A * beta - y_meas; %calculates residuals
    SSQ = r'*r; % sum of squared residuals
    SStotal = (length(y_meas)-1) * var(y_meas);
    rsq = 1 - SSQ/SStotal;
    text(x_fit(1),y_fit(1),['R^2 = ' num2str(rsq)],'FontSize',12,'Color',autumn(1),...
        'VerticalAlignment','Bottom')
    
    % Linear regression (for glycogen EXCESS)
    x = Gly_Cfed_exp(Gly_Cfed_exp>0); y_meas = AcCoA_Cfed_exp(Gly_Cfed_exp>0);
    A = [x ones(size(x))];
    s_y = AcCoA_Cfed_exp_STD(Gly_Cfed_exp>0);
    S_y = diag( s_y.^2 ); % covariance matrix
    beta = inv( A'*inv(S_y)*A) * A' *inv(S_y) * y_meas;
    
    x_fit = 0:0.1:0.4;
    y_fit = [x_fit' ones(size(x_fit'))] * beta;
    plot(x_fit,y_fit,'Color',summer(1));
    
    r = A * beta - y_meas; %calculates residuals
    SSQ = r'*r; % sum of squared residuals
    SStotal = (length(y_meas)-1) * var(y_meas);
    rsq = 1 - SSQ/SStotal;
    text(x_fit(end),y_fit(end),['R^2 = ' num2str(rsq)],'FontSize',12,'Color',summer(1),...
        'VerticalAlignment','Bottom','HorizontalAlignment','Right')
end

% pretty plot
lgd = legend([p1 p2 p3 p4 p5 p6],'Unidentified','PAO','PAO I', 'PAO II', 'GAO','EBPR sludge',...
    'Location','East'); legend('boxoff');
lgd.FontSize = 14;
title(lgd,'Enrichment type')
lgd.Position = [0.75    0.38    0.2214    0.3107];
set(gca,'color','white')

yl_main=ylim; xl_main=xlim;

% xticks([(0-6/18)	(3/15-6/18)     0	(9/21-6/18)     (12/24-6/18)    (18/30-6/18)])
% xticklabels({'0Gly:6Ac','0.5Gly:6Ac','1Gly:6Ac','1.5Gly:6Ac','2Gly:6Ac','3Gly:6Ac'})

xticks([(0-6/18)	(1/6-6/18)     0	(1/2-6/18)     (2/3-6/18)    (5/6-6/18)])
xticklabels({'0 Gly:1 Ac','^{1}/_{6}Gly:^{5}/_{6}Ac','^{2}/_{6}Gly:^{4}/_{6}Ac','^{3}/_{6}Gly:^{3}/_{6}Ac','^{4}/_{6}Gly:^{2}/_{6}Ac','^{5}/_{6}Gly:^{1}/_{6}Ac'})
ax.FontSize = 12;

% make insert of whole triangle
% create smaller axes in top right, and plot on it
axes('Position',[.705 .7 .2 .2])

p1 = plot(Gly_Cfed,theoreticalY_PHB_Cfed*ones(size(Gly_Cfed)),':k', 'linewidth', 2,'Color',[0.75 0.75 0.75]);
hold on
p2 = plot(Gly_Cfed,AcCoA_Cfed_optimum,'--k', 'linewidth', 1);
p3 = plot(Gly_Cfed,AcCoA_Cfed_minimum, ':k', 'linewidth', 2);
p4 = fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
    [AcCoA_Cfed_minimum;AcCoA_Cfed_optimum(end:-1:1);AcCoA_Cfed_minimum(1)],...
    'y','FaceColor',cmap(end-5,:),'FaceAlpha',0.2,'LineStyle','none');

% get rid of standard axes decorations
set(gca, 'box', 'off')

% Fix the axes sizes
y_axis_size = 1.04;
% x_axis_size = 0.8;
axis([-0.5+0.01 x_axis_size-0.01 0 y_axis_size])
yl=ylim; xl=xlim;
ax = gca;
ax.YAxisLocation = 'origin';
ax.ColorOrderIndex = 1;

hold on
% PLOT experimental data
p = plot(Gly_Cfed_exp,AcCoA_Cfed_exp,'s','Color',[0.5 0.5 0.5],'MarkerSize',3);
p = plot(Gly_Cfed_exp(PAO_idx),AcCoA_Cfed_exp(PAO_idx),'o','MarkerSize',3);
p = plot(Gly_Cfed_exp(PAO1_idx),AcCoA_Cfed_exp(PAO1_idx),'.','MarkerSize',6);
p = plot(Gly_Cfed_exp(PAO2_idx),AcCoA_Cfed_exp(PAO2_idx),'.','MarkerSize',3);
p = plot(Gly_Cfed_exp(GAO_idx),AcCoA_Cfed_exp(GAO_idx),'xk','linewidth',0.5,'MarkerSize',6);
p = plot(Gly_Cfed_exp(EBPR_idx),AcCoA_Cfed_exp(EBPR_idx),'.','MarkerSize',6);

set(gca,'color','white')

p = fill([xl_main(1) xl_main(1) xl_main(2) xl_main(2)],...
    [yl_main yl_main(end:-1:1)],...
    'k','FaceAlpha',0.2,'LineStyle','none');

lgd = legend([p1 p2 p3 p4],'Theoretical max','Optimum','Minimum','Solution space',...
    'Location','NorthEast'); legend('boxoff');
lgd.FontSize = 10;
lgd.Position = [0.865-0.001 0.795+0.026 0.0764 0.0427];

xticks([(0-6/18)	0     (2/3-6/18)    (1-6/18)])
xticklabels({'0 Gly:1 Ac','^{2}/_{6}Gly:^{4}/_{6}Ac','^{4}/_{6}Gly:^{2}/_{6}Ac','1 Gly:0 Ac'})

ax.TickLength = [0.02 0.035];

%% PLOT PrCoA YIELD
figure('Name','PrCoA','NumberTitle','off'), clf, hold on

theoreticalY_PH2MV_Cfed = 4/5; %cmolPH2MV/cmolAc
plot(Gly_Cfed,theoreticalY_PH2MV_Cfed*ones(size(Gly_Cfed)),':k', 'linewidth', 2,'Color',[0.75 0.75 0.75])

plot(Gly_Cfed(Gly_Cfed<=0),PrCoA_Cfed_optimum(Gly_Cfed<=0),'--r','Color',[226/256 26/256 26/256], 'linewidth', 2)
plot(Gly_Cfed(Gly_Cfed>=0),PrCoA_Cfed_optimum(Gly_Cfed>=0),'--g','Color',[165/256 202/256 26/256], 'linewidth', 2)
plot(Gly_Cfed,PrCoA_Cfed_minimum, ':k', 'linewidth', 2)

fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
    [PrCoA_Cfed_minimum;PrCoA_Cfed_optimum(end:-1:1);PrCoA_Cfed_minimum(1)],...
    'y','FaceColor',cmap(end-5,:),'FaceAlpha',0.1,'LineStyle','none')

plot(0,0,'p','Color', [0/256 166/256 214/256],'MarkerSize',14, 'linewidth', 2)

if Yagci_model
    plot(Gly_Cfed_Yagci_PAO,PrCoA_Cfed_Yagci_PAO,':r','Color',[109/256 23/256 127/256], 'linewidth', 2)
    plot(Gly_Cfed_Yagci_GAO,PrCoA_Cfed_Yagci_GAO,':g','Color',[0/256 136/256 145/256], 'linewidth', 2)
    
    text(-0.268, 0.375, {'Yagci (2003) PAO model'}, 'fontsize', 14,'rotation', -56,...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom','Color',[109/256 23/256 127/256])
    text(0.37, 0.26, {'Yagci (2003) GAO model'}, 'fontsize', 14,'rotation', +36,...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom','Color', [0/256 136/256 145/256])
end

if H2_TCA
    H2_color = [0/256 136/256 145/256];
    plot(Gly_Cfed,PrCoA_Cfed_optimum_H2,'--r','Color',H2_color, 'linewidth', 2)
    plot(Gly_Cfed,PrCoA_Cfed_minimum_H2,':g','Color',H2_color, 'linewidth', 2)
    fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
        [PrCoA_Cfed_minimum_H2;PrCoA_Cfed_optimum_H2(end:-1:1);PrCoA_Cfed_minimum_H2(1)],...
        'y','FaceColor',H2_color,'FaceAlpha',0.1,'LineStyle','none');
    text(0.27, 0.1, {'with','H_2', 'prod.'}, 'fontsize', 14,'rotation', -0,...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom','Color',H2_color)
    
    TCA_color = [109/256 23/256 127/256];
    plot(Gly_Cfed,PrCoA_Cfed_optimum_TCA,'--r','Color',TCA_color, 'linewidth', 2)
    plot(Gly_Cfed,PrCoA_Cfed_minimum_TCA,':g','Color',TCA_color, 'linewidth', 2)
    fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
        [PrCoA_Cfed_minimum_TCA;PrCoA_Cfed_optimum_TCA(end:-1:1);PrCoA_Cfed_minimum_TCA(1)],...
        'y','FaceColor',TCA_color,'FaceAlpha',0.1,'LineStyle','none');
    
    text(-0.3, 0.2, {'with full TCA'}, 'fontsize', 14,'rotation', -0,...
        'HorizontalAlignment','left', 'VerticalAlignment','bottom','Color',TCA_color)
end

% get rid of standard axes decorations
set(gca, 'box', 'off')

% Fix the axes sizes
y_axis_size = 0.43;
% x_axis_size = 0.8;
axis([-0.5+0.11 x_axis_size-0.41 0 y_axis_size])
yl=ylim; xl=xlim;
ax = gca;
ax.YAxisLocation = 'origin';
ax.ColorOrderIndex = 1;

% the static texts
text(0,    y_axis_size, {'PrCoA* accumulated per C consumed','(Cmol/Cmol)'}, 'rotation', 0, 'fontsize', 14,...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom')

text(xl(2) , yl(1)-0.02, 'Stoichiometric proportion C sources consumed (Cmol)', 'fontsize', 14,...
    'HorizontalAlignment','right', 'VerticalAlignment','top')

text(0.36, yl(2)-0.158, {'C) Glycolysis +','redTCA branch'}, 'fontsize', 14, 'rotation', +36, ...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','Color', [165/256 202/256 26/256])

text(-0.25, yl(2)-0.042, {'A) Glyoxylate Shunt +','Glycolysis'}, 'fontsize', 14, 'rotation', -56, ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom','Color',[226/256 26/256 26/256])

text(0.01, 0.01, {'B) Glycolysis only'}, 'fontsize', 14,'rotation', +0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [0/256 166/256 214/256])

% PLOT experimental data
if errorbars
    errorbar(Gly_Cfed_exp,PrCoA_Cfed_exp,PrCoA_Cfed_exp_STD,PrCoA_Cfed_exp_STD,Gly_Cfed_exp_STD,Gly_Cfed_exp_STD,...
        'Color',[0.8 0.8 0.8],'LineStyle','none','LineWidth',0.2)
end

leg = {strcat(num2str(1), '.  ', Author{1},' (',Year{1}, ')')};
counter = 1;
pos=1;
for ii = 2:length(Author)
    if isequal(Author{ii},Author{ii-1}) && isequal(Year{ii},Year{ii-1})
        leg{pos} = strcat(num2str(counter), '-',num2str(ii), '.  ', Author{ii},' (',Year{ii}, ')');
    else
        pos = pos+1;
        counter = ii;
        leg{pos} = strcat(num2str(ii), '.  ', Author{ii},' (',Year{ii}, ')');
    end
end


enrichment = TXT(2:end,3);
PAO_idx = contains(enrichment,'PAO');
PAO1_idx = contains(enrichment,'PAO I ');
PAO2_idx = contains(enrichment,'PAO II');
GAO_idx = contains(enrichment,'GAO');
EBPR_idx = contains(enrichment,'EBPR');
%unidentified = contains(enrichment,'????');


p1 = plot(Gly_Cfed_exp,PrCoA_Cfed_exp,'s','Color',[0.5 0.5 0.5],'MarkerSize',6);
p2 = plot(Gly_Cfed_exp(PAO_idx),PrCoA_Cfed_exp(PAO_idx),'o','MarkerSize',6);
p3 = plot(Gly_Cfed_exp(PAO1_idx),PrCoA_Cfed_exp(PAO1_idx),'.','MarkerSize',12);
p4 = plot(Gly_Cfed_exp(PAO2_idx),PrCoA_Cfed_exp(PAO2_idx),'.','MarkerSize',6);
p5 = plot(Gly_Cfed_exp(GAO_idx),PrCoA_Cfed_exp(GAO_idx),'xk','linewidth',1,'MarkerSize',12);
p6 = plot(Gly_Cfed_exp(EBPR_idx),PrCoA_Cfed_exp(EBPR_idx),'.','MarkerSize',12);


text(Gly_Cfed_exp,PrCoA_Cfed_exp,num2str([1:length(PrCoA_Cfed_exp)]'),...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom')

text( xl(1)+0.04 , y_axis_size , leg , ...
    'HorizontalAlignment','right', 'VerticalAlignment','top',...
    'FontSize',12)

if regression
    % Linear regression (for glycogen LIMITATION)
    x = Gly_Cfed_exp(Gly_Cfed_exp<0); y_meas = PrCoA_Cfed_exp(Gly_Cfed_exp<0);
    A = [x ones(size(x))];
    s_y = PrCoA_Cfed_exp_STD(Gly_Cfed_exp<0);
    S_y = diag( s_y.^2 ); % covariance matrix
    beta = inv( A'*inv(S_y)*A) * A' *inv(S_y) * y_meas;
    
    x_fit = -0.4:0.1:0;
    y_fit = [x_fit' ones(size(x_fit'))] * beta;
    plot(x_fit,y_fit,'Color',autumn(1));
    
    r = A * beta - y_meas; %calculates residuals
    SSQ = r'*r; % sum of squared residuals
    SStotal = (length(y_meas)-1) * var(y_meas);
    rsq = 1 - SSQ/SStotal;
    text(x_fit(1),y_fit(1),['R^2 = ' num2str(rsq)],'FontSize',12,'Color',autumn(1),...
        'VerticalAlignment','Bottom')
    
    % Linear regression (for glycogen EXCESS)
    x = Gly_Cfed_exp(Gly_Cfed_exp>0); y_meas = PrCoA_Cfed_exp(Gly_Cfed_exp>0);
    A = [x ones(size(x))];
    s_y = PrCoA_Cfed_exp_STD(Gly_Cfed_exp>0);
    S_y = diag( s_y.^2 ); % covariance matrix
    beta = inv( A'*inv(S_y)*A) * A' *inv(S_y) * y_meas;
    
    x_fit = 0:0.1:0.4;
    y_fit = [x_fit' ones(size(x_fit'))] * beta;
    plot(x_fit,y_fit,'Color',summer(1));
    
    r = A * beta - y_meas; %calculates residuals
    SSQ = r'*r; % sum of squared residuals
    SStotal = (length(y_meas)-1) * var(y_meas);
    rsq = 1 - SSQ/SStotal;
    text(x_fit(end),y_fit(end),['R^2 = ' num2str(rsq)],'FontSize',12,'Color',summer(1),...
        'VerticalAlignment','Bottom','HorizontalAlignment','Right')
end

% pretty plot
lgd = legend([p1 p2 p3 p4 p5 p6],'Unidentified','PAO','PAO I', 'PAO II', 'GAO','EBPR sludge',...
    'Location','SouthEast'); legend('boxoff');
lgd.FontSize = 14;
title(lgd,'Enrichment type')
lgd.Position = [0.75 0.2 0.2214 0.3107];
set(gca,'color','white')

yl_main=ylim; xl_main=xlim;

xticks([(0-6/18)	(1/6-6/18)     0	(1/2-6/18)     (2/3-6/18)    (5/6-6/18)])
xticklabels({'0 Gly:1 Ac','^{1}/_{6}Gly:^{5}/_{6}Ac','^{2}/_{6}Gly:^{4}/_{6}Ac','^{3}/_{6}Gly:^{3}/_{6}Ac','^{4}/_{6}Gly:^{2}/_{6}Ac','^{5}/_{6}Gly:^{1}/_{6}Ac'})
ax.FontSize = 12;

% make insert of whole triangle
% create smaller axes in top right, and plot on it
axes('Position',[.705 .75 .2 .2])

p1 = plot(Gly_Cfed,theoreticalY_PHB_Cfed*ones(size(Gly_Cfed)),':k', 'linewidth', 2,'Color',[0.75 0.75 0.75]);
hold on
p2 = plot(Gly_Cfed,PrCoA_Cfed_optimum,'--k', 'linewidth', 1);
p3 = plot(Gly_Cfed,PrCoA_Cfed_minimum, ':k', 'linewidth', 2);
p4 = fill([Gly_Cfed;Gly_Cfed(end:-1:1);Gly_Cfed(1)],...
    [PrCoA_Cfed_minimum;PrCoA_Cfed_optimum(end:-1:1);PrCoA_Cfed_minimum(1)],...
    'y','FaceColor',cmap(end-5,:),'FaceAlpha',0.2,'LineStyle','none');

% get rid of standard axes decorations
set(gca, 'box', 'off')

% Fix the axes sizes
y_axis_size = 1.04;
% x_axis_size = 0.8;
axis([-0.5+0.01 x_axis_size-0.01 0 y_axis_size])
yl=ylim; xl=xlim;
ax = gca;
ax.YAxisLocation = 'origin';
ax.ColorOrderIndex = 1;

hold on
% PLOT experimental data
p = plot(Gly_Cfed_exp,PrCoA_Cfed_exp,'s','Color',[0.5 0.5 0.5],'MarkerSize',3);
p = plot(Gly_Cfed_exp(PAO_idx),PrCoA_Cfed_exp(PAO_idx),'o','MarkerSize',3);
p = plot(Gly_Cfed_exp(PAO1_idx),PrCoA_Cfed_exp(PAO1_idx),'.','MarkerSize',6);
p = plot(Gly_Cfed_exp(PAO2_idx),PrCoA_Cfed_exp(PAO2_idx),'.','MarkerSize',3);
p = plot(Gly_Cfed_exp(GAO_idx),PrCoA_Cfed_exp(GAO_idx),'xk','linewidth',0.5,'MarkerSize',6);
p = plot(Gly_Cfed_exp(EBPR_idx),PrCoA_Cfed_exp(EBPR_idx),'.','MarkerSize',6);

set(gca,'color','white')

p = fill([xl_main(1) xl_main(1) xl_main(2) xl_main(2)],...
    [yl_main yl_main(end:-1:1)],...
    'k','FaceAlpha',0.2,'LineStyle','none');

lgd = legend([p1 p2 p3 p4],'Theoretical max','Optimum','Minimum','Solution space',...
    'Location','NorthEast'); legend('boxoff');
lgd.FontSize = 10;
lgd.Position = [0.865-0.001 0.795-0.026 0.0764 0.0427];

% xticks([(0-6/18)	0	(18/30-6/18)    (72/84-6/18)])
% xticklabels({'0Gly:6Ac','1Gly:6Ac','3Gly:6Ac','12Gly:6Ac'})

xticks([(0-6/18)	0     (2/3-6/18)    (1-6/18)])
xticklabels({'0 Gly:1 Ac','^{2}/_{6}Gly:^{4}/_{6}Ac','^{4}/_{6}Gly:^{2}/_{6}Ac','1 Gly:0 Ac'})

ax.TickLength = [0.02 0.035];

%% ------------- FLUX through PATHWAYS -------------
%OPTIMUM SCENARIO
figure('Name','Optimum scenario','NumberTitle','off'), clf, hold on

sum_pathways = Glycolysis_Cfed_optimum + TCAGOX_Cfed_optimum + oxTCA_Cfed_optimum + PDH_Cfed_optimum;
plot(Gly_Cfed,sum_pathways,':m', 'linewidth', 2,'Color',[0.75 0.75 0.75])
sum_sink_pathways = -(1/2)*AcCoA_Cfed_optimum -(1/3)*PrCoA_Cfed_optimum + redTCA_Cfed_optimum;
plot(Gly_Cfed,sum_sink_pathways,':k', 'linewidth', 2,'Color',[0.75 0.75 0.75])

%plot(Gly_Cfed,Glycolysis_Cfed_optimum ,'--r','Color', [0/256 166/256 214/256],'MarkerSize',14, 'linewidth', 2)

y = [Glycolysis_Cfed_optimum';...
    TCAGOX_Cfed_optimum';...
    oxTCA_Cfed_optimum';...
    PDH_Cfed_optimum'];
h = area(Gly_Cfed,y','LineStyle','none');
h(1).FaceColor = [0/256 166/256 214/256];
h(1).FaceAlpha = 0.2;
h(2).FaceColor = [226/256 26/256 26/256];
h(2).FaceAlpha = 0.2;
h(3).FaceColor = [225/256 196/256 0/256];
h(3).FaceAlpha = 0.2;
h(4).FaceColor = [107/256 134/256 137/256];
h(4).FaceAlpha = 0.2;

z = [redTCA_Cfed_optimum';...
    -(1/2)*AcCoA_Cfed_optimum';...
    -(1/3)*PrCoA_Cfed_optimum'];
h = area(Gly_Cfed,z','LineStyle','none');
h(1).EdgeColor = [165/256 202/256 26/256];
h(1).FaceColor = [165/256 202/256 26/256];
h(1).FaceAlpha = 0.2;
h(2).FaceColor = [109/256 23/256 127/256];
h(2).FaceAlpha = 0.2;
h(3).FaceColor = [230/256 70/256 22/256];
h(3).FaceAlpha = 0.2;

% get rid of standard axes decorations
set(gca, 'box', 'off')

% Fix the axes sizes
y_axis_size = 1.04;
% x_axis_size = 0.8;
axis([-0.5+0.11 x_axis_size-0.41 -y_axis_size y_axis_size])
yl=ylim; xl=xlim;
ax = gca;
ax.YAxisLocation = 'origin';
ax.XAxisLocation = 'origin';
ax.ColorOrderIndex = 1;

% Flux modes
errorbar(-1/6, 0.9, 0.15, 'horizontal','Color',[226/256 26/256 26/256], 'linewidth', 2)
text(-1/6, 0.9, 'A)','Color',[226/256 26/256 26/256], 'fontsize', 14,...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom')
errorbar(1/6, 0.9, 0.15, 'horizontal','Color',[165/256 202/256 26/256], 'linewidth', 2)
text(1/6, 0.9, 'C)','Color',[165/256 202/256 26/256], 'fontsize', 14,...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom')
plot(0,0.9,'p','Color', [0/256 166/256 214/256],'MarkerSize',14, 'linewidth', 2)
text(0.01 ,0.9, 'B)','Color',[0/256 166/256 214/256], 'fontsize', 14,...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom')

% the static texts
text(0,    y_axis_size, {'Electrons harvested/sinked per C consumed','(Emol/Cmol)'}, 'rotation', 0, 'fontsize', 14,...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom')

text(xl(1)-0.05 , 0, {'Stoichiometric proportion','C sources consumed','(Cmol)'}, 'fontsize', 14,...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', 'BackgroundColor', [1 1 1])

text(xl(2)-0.02, -0.25, {'redTCA branch'}, 'fontsize', 14, 'rotation', 0, ...
    'HorizontalAlignment','right', 'VerticalAlignment','bottom','Color', [165/256 202/256 26/256])

text(xl(1)+0.1, 0.2, {'Glyoxylate Shunt'}, 'fontsize', 14, 'rotation', -0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color',[226/256 26/256 26/256])

text(0.2, 0.2, {'Glycolysis'}, 'fontsize', 14,'rotation', +0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [0/256 166/256 214/256])

text(0.2, 0.5, {'PyrDH'}, 'fontsize', 14,'rotation', +0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [107/256 134/256 137/256])

text(0.1, -0.35, {'AcCoA* acc.'}, 'fontsize', 14,'rotation', +0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [109/256 23/256 127/256])

text(-0.31, -0.28, {'PrCoA* acc.'}, 'fontsize', 14,'rotation', +0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [230/256 70/256 22/256])

text(xl(2), 0.82, {'Sum electrons harvested'}, 'fontsize', 14,'rotation', +10, ...
    'HorizontalAlignment','right', 'VerticalAlignment','middle','Color', [0.75 0.75 0.75])

text(xl(2), -0.82, {'Sum electrons sinked'}, 'fontsize', 14,'rotation', -10, ...
    'HorizontalAlignment','right', 'VerticalAlignment','middle','Color', [0.75 0.75 0.75])

xticks([(0-6/18)	(1/6-6/18)     0	(1/2-6/18)     (2/3-6/18)    (5/6-6/18)])
xticklabels({'0 Gly:1 Ac','^{1}/_{6}Gly:^{5}/_{6}Ac','^{2}/_{6}Gly:^{4}/_{6}Ac','^{3}/_{6}Gly:^{3}/_{6}Ac','^{4}/_{6}Gly:^{2}/_{6}Ac','^{5}/_{6}Gly:^{1}/_{6}Ac'})
ax.FontSize = 12;

%MINIMUM SCENARIO
figure('Name','Minimum scenario','NumberTitle','off'), clf, hold on

sum_pathways = Glycolysis_Cfed_minimum + TCAGOX_Cfed_minimum + oxTCA_Cfed_minimum + PDH_Cfed_minimum;
plot(Gly_Cfed,sum_pathways,':m', 'linewidth', 2,'Color',[0.75 0.75 0.75])
sum_sink_pathways = -(1/2)*AcCoA_Cfed_minimum -(1/3)*PrCoA_Cfed_minimum + redTCA_Cfed_minimum;
plot(Gly_Cfed,sum_sink_pathways,':k', 'linewidth', 2,'Color',[0.75 0.75 0.75])

%plot(Gly_Cfed,Glycolysis_Cfed_minimum ,'--r','Color', [0/256 166/256 214/256],'MarkerSize',14, 'linewidth', 2)

y = [Glycolysis_Cfed_minimum';...
    TCAGOX_Cfed_minimum';...
    oxTCA_Cfed_minimum';...
    PDH_Cfed_minimum'];
h = area(Gly_Cfed,y','LineStyle','none');
h(1).FaceColor = [0/256 166/256 214/256];
h(1).FaceAlpha = 0.2;
h(2).FaceColor = [226/256 26/256 26/256];
h(2).FaceAlpha = 0.2;
h(3).FaceColor = [225/256 196/256 0/256];
h(3).FaceAlpha = 0.2;
h(4).FaceColor = [107/256 134/256 137/256];
h(4).FaceAlpha = 0.2;

z = [redTCA_Cfed_minimum';...
    -(1/2)*AcCoA_Cfed_minimum';...
    -(1/3)*PrCoA_Cfed_minimum'];
h = area(Gly_Cfed,z','LineStyle','none');
h(1).EdgeColor = [165/256 202/256 26/256];
h(1).FaceColor = [165/256 202/256 26/256];
h(1).FaceAlpha = 0.2;
h(2).FaceColor = [109/256 23/256 127/256];
h(2).FaceAlpha = 0.2;
h(3).FaceColor = [230/256 70/256 22/256];
h(3).FaceAlpha = 0.2;

% get rid of standard axes decorations
set(gca, 'box', 'off')

% Fix the axes sizes
y_axis_size = 1.04;
% x_axis_size = 0.8;
axis([-0.5+0.11 x_axis_size-0.41 -y_axis_size y_axis_size])
yl=ylim; xl=xlim;
ax = gca;
ax.YAxisLocation = 'origin';
ax.XAxisLocation = 'origin';
ax.ColorOrderIndex = 1;


% the static texts
text(0,    y_axis_size, {'Electrons harvested/sinked per C consumed','(Emol/Cmol)'}, 'rotation', 0, 'fontsize', 14,...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom')

text(xl(1)-0.05 , 0, {'Stoichiometric proportion','C sources consumed','(Cmol)'}, 'fontsize', 14,...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', 'BackgroundColor', [1 1 1])

text(0.2, -0.3, {'redTCA branch'}, 'fontsize', 14, 'rotation', 0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','bottom','Color', [165/256 202/256 26/256])

text(-0.2, 0.3, {'Glyoxylate Shunt'}, 'fontsize', 14, 'rotation', -0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color',[226/256 26/256 26/256])

text(0.2, 0.2, {'Glycolysis'}, 'fontsize', 14,'rotation', +0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [0/256 166/256 214/256])

text(0.2, 0.7, {'PyrDH'}, 'fontsize', 14,'rotation', +0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [107/256 134/256 137/256])

text(-0.31, -0.13, {'AcCoA* acc.'}, 'fontsize', 14,'rotation', +0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [109/256 23/256 127/256])

text(-0.31, -0.28, {'PrCoA* acc.'}, 'fontsize', 14,'rotation', +0, ...
    'HorizontalAlignment','left', 'VerticalAlignment','middle','Color', [230/256 70/256 22/256])

text(xl(2), 0.91, {'Sum electrons harvested'}, 'fontsize', 14,'rotation', +7, ...
    'HorizontalAlignment','right', 'VerticalAlignment','middle','Color', [0.75 0.75 0.75])

text(xl(2), -0.91, {'Sum electrons sinked'}, 'fontsize', 14,'rotation', -7, ...
    'HorizontalAlignment','right', 'VerticalAlignment','middle','Color', [0.75 0.75 0.75])

xticks([(0-6/18)	(1/6-6/18)     0	(1/2-6/18)     (2/3-6/18)    (5/6-6/18)])
xticklabels({'0 Gly:1 Ac','^{1}/_{6}Gly:^{5}/_{6}Ac','^{2}/_{6}Gly:^{4}/_{6}Ac','^{3}/_{6}Gly:^{3}/_{6}Ac','^{4}/_{6}Gly:^{2}/_{6}Ac','^{5}/_{6}Gly:^{1}/_{6}Ac'})
ax.FontSize = 12;

