%% FBA
% Written by: 
% Leonor Guedes da Silva
% LeonorGuedesdaSilva@gmail.com
% 
% Runs in MATLAB R2018a for macOSX
% 
% Package requirements:
% - xlwrite (alternative available in Windows systems: xlswrite)
% https://nl.mathworks.com/matlabcentral/fileexchange/38591-xlwrite-generate-xls-x-files-without-excel-on-mac-linux-win
% 
% Input:
%       - Stoichiometrix matrix (in mol)
% 
% Output:
%       - Flux distributions (in mol/Cmol consumed)
%
%%
clearvars
close all
clc

%% Import stoichiometric matrix from Excel
T = readtable('Stoichiometry_FBA.csv');

N = T{1:end-1,2:end};
rev = T{end,2:end};
lab_met = T{1:end-1,1}';
lab_flx = T.Properties.VariableNames(2:end);

%Generate variables with metabolite and flux names associated with their
%position in the stoichiometric matrix
for ii = 1:length(lab_met)
   v = genvarname(lab_met{ii});
   eval([v ' = ii']);
end
for ii = 1:length(lab_flx)
   v = genvarname(lab_flx{ii}); %#ok<*DEPGENAM>
   eval([v ' = ii']);
end

%% RESULTS FILENAME

results_filename = 'S_results_FBA.xls';


%% ------------ Simulations ------------
sim = 'AcCoAPrCoAprod';
max_or_min = 'min'; %maximum simulation: "max"; minimum simulation: "min"
target = CO2Prod; %flux to max or minimize

% DEFAULT SETTINGS
meas_idx = [GlycDeg AcUpt];
upperbounded_idx = [];
upperbounded_flx = [ ]'; %column vector
knockouts = [PHBAcc PHVAcc NOG SucDH TCA H2Prod PPP];  %sets KO's to ZERO = OFF

    
switch sim
    case 'AcCoAPrCoAprod'
        disp('------------Optimization for anaerobic PHA from acetate and glycogen------------')
        disp('------------(without full TCA nor H2 production)------------')
        
        knockouts = [PHBAcc PHVAcc NOG SucDH TCA H2Prod PPP];

    case 'fullTCA'
        disp('------------Optimization for anaerobic PHA from acetate and glycogen------------')
        disp('------------(with full TCA)------------')
       
        knockouts = [PHBAcc PHVAcc NOG SucDH H2Prod PPP];

    case 'H2Prod'
        disp('------------Optimization for anaerobic PHA from acetate and glycogen------------')
        disp('------------(with H2 production)------------')

        knockouts = [PHBAcc PHVAcc NOG TCA SucDH PPP];
    
    %% ------------ EXTRA Simulations ------------
    case 'PHB_PHV'
        disp('------------Optimization for anaerobic PHB+PHV from acetate and glycogen------------')
        disp('------------(without full TCA nor H2 production)------------')

        knockouts = [AcCoAProd PrCoAProd NOG SucDH TCA H2Prod PPP];

    case 'limitedPEPC_H2Prod'
        disp('------------What if PAOs are kinetically limited in their PEPC+redTCA route?------------')
        disp('------------(with H2 production allowed)------------') 
        
        knockouts = [PHBAcc PHVAcc NOG SucDH TCA PPP];
        upperbounded_idx = PEPC;
        upperbounded_flx = 0.05;
        
    case 'succinate_export'
        disp('------------What if cells excrete succinate?------------')
        disp('------------(without full TCA nor H2 production)------------')
        
        % add succinate exp reaction
        lab_flx{end+1} = 'SuccExp';
        N(SucCoA,end+1) = -1;
        rev(end+1) = 0;
        
        knockouts = [AcCoAProd PrCoAProd NOG SucDH TCA H2Prod PPP];
        %target = PHBAcc;
        
    case 'H2supplied'
        disp('------------What if cells get H2 from methanogenes?------------')
        disp('------------(without full TCA but with H2 consumption)------------')
        
        % change direction of H2Prod
        N(:,H2Prod) = -1 * N(:,H2Prod);
        
        knockouts = [AcCoAProd PrCoAProd NOG SucDH TCA PPP];
        % knockouts = [PHBAcc PHVAcc NOG SucDH TCA PPP];
        target = PHBAcc;

    case 'GlycineProd'
        disp('------------What if cells also accumulate glycine?------------')
        disp('------------(without full TCA nor H2 production)------------')
        
        % add glycine exp reaction
        lab_flx{end+1} = 'GlycineAcc';
        N(PYR,end+1) = -1;
        rev(end+1) = 0;
        rev(PDC) = 1; %make PDH reversible = PFOR
        
        knockouts = [AcCoAProd PrCoAProd NOG SucDH TCA H2Prod PPP];
        
    case 'GlutamateProd'
        disp('------------What if cells also accumulate glutamate (+ NH3 consumption)?------------')
        disp('------------(without full TCA nor H2 production)------------')
        
        % add glutamate exp reaction
        lab_flx{end+1} = 'GlutamateAcc'; % aKG + NH3 + NADPH -> Glutamate
        N(SucCoA,end+1) = -1;
        N(electrons,end) = -2 - 2; %consumes NADPH! and SucCoA to aKG also
        N(CO2,end) = -1; %consumes 1 CO2 from SucCoA to aKG
        rev(end+1) = 0;
        rev(PDC) = 1; %make PDH reversible = PFOR
        
        knockouts = [AcCoAProd PrCoAProd NOG SucDH TCA H2Prod PPP];
        
    case 'NoAcetate'
        disp('------------What if.. no acetate?------------')
        disp('------------(with full TCA and H2 production)------------')
        
        knockouts = [PHBAcc PHVAcc NOG SucDH PPP AcUpt];
        
    otherwise
        disp('select a valid simulation scenario')
end


%% Optimization specs
if strcmp(max_or_min, 'max'), sign = 1; 
elseif strcmp(max_or_min, 'min'), sign = -1; end

objective = sign * target; 

disp(['Turned-off reactions: ' strjoin( lab_flx(knockouts) , ', ')])

%% Network analysis
network = setdiff(1:size(N,2),knockouts);

transportR =  [meas_idx target]; % transport reactions in N
[~, transportR_network] = intersect(lab_flx(network),lab_flx(transportR));

network_analysis( N(:,network), lab_flx(network) , lab_met, rev(network), transportR_network) 

%% Write results to Excel
sheetname = strcat(sim, '_', max_or_min, lab_flx{target} );
disp(['Sheetname: ' sheetname])
xlwrite(results_filename,[{'Objective','Knockouts','f'},lab_flx], sheetname, 'A1')
row_no = 2;

%% Loop over different proportions of acetate:glycogen consumed
for f = [0 1/6 2/6 3/6 4/6 5/6 1] %Cmol acetate / Cmol consumed
    
    Acupt_flx = f; %Cmol acetate / Cmol consumed
    Glyc = 1 - f; %Cmol glycogen / Cmol consumed
    
    Acupt_flx_mol = Acupt_flx /2; % mol acetate / Cmol consumed
    Glyc_mol = Glyc / 6; % mol glycogen / Cmol consumed
    
    meas_flx = [Glyc_mol; Acupt_flx_mol]; %column vector    
    
    % Perform optimization
    v_opt = FBAopt(N,objective,meas_idx,meas_flx,...
                   knockouts,rev,upperbounded_idx,upperbounded_flx);
   
    % write results in excel
    text = {lab_flx(abs(objective)), strjoin( lab_flx(knockouts) , ', ')};  
    xlwrite(results_filename,text, sheetname, strcat('A',num2str(row_no)))

    num = [ f, v_opt'];
    xlwrite(results_filename, num, sheetname, strcat('C',num2str(row_no))); 
    row_no = row_no + 1;
    
end
