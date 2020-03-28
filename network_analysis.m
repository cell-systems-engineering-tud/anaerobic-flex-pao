function network_analysis( N, lab_flx, lab_met, rev, transportR)  

%% Print reactions
print_reactions_rev( N, lab_flx, lab_met, rev )     

%% Model specs:
disp('------------Model specs:------------')
numberof_transportR = length(transportR);
disp( [ '# Transport reactions (In-/Outputs): ' num2str(numberof_transportR)] );
intraR =  size(N,2)-numberof_transportR;
disp( [ '# Intracellular reactions: ' num2str(intraR)] );
metbal =  size(N,1);
disp( [ '# Metabolite balances: ' num2str(metbal)] );
dof =  size(N,2)- rank(N);
disp( [ 'Degrees of freedom: ' num2str(dof)] );
CMs = size(N,1) - rank(N);
disp( [ '# Conserved moieties: ' num2str(CMs)] );

%------------Which are the conserved moieties?------------
if CMs~=0
    NS = null( N', 'r' );
    for ii=1:size(NS,2)
        CM_idx = find(NS(:,ii));
        disp(strjoin( lab_met(CM_idx) , ' + '));
    end
end

%------------Blocked reactions?------------
NS = null( N, 'r' );
blocked_reactions = [];
for ii=1:length( NS ),
    if all( abs( NS(ii,:) ) < 1e-12 ), % entries < 1e-12?
        %fprintf( 'Reaction %s is blocked\n', lab_flx{ ii } );
        blocked_reactions(end+1) = ii;
    end
end
if isempty(blocked_reactions)
    disp( 'Blocked reactions?: No' )
else disp([ 'Blocked reactions?: ' strjoin( lab_flx(blocked_reactions) , ', ') ])
end

%%
% Test the model for intracellular cycles ? are there reactions that can
% form cycles without an influence on the extracellular measurements?
disp('------------Internal cycles?------------')
%Constrain measured fluxes
disp('Measured fluxes:')
disp(strjoin(lab_flx(transportR) , ', '))
C = eye(size(N,2)); %same size as the number of fluxes
C = C(transportR,:); %adds a constraint in the measured fluxes (all transport fluxes)

NS = round(null( [N;C], 'r' ), 12); %round to prevent innacuracies
% disp('Nullspace:')
% disp(NS)
disp(['Internal cycles found: ' num2str(size(NS,2))])

invalid_internal_cycles = 0; %counter
% are there fluxes in the nullspace?
disp( 'Fluxes that cannot be determined:' );
for ii=1:size( NS, 2 )
    % is the internal cycle valid (i.e. does it respect irreversibilities?)
    NS_irrev = NS(:,ii) ./ abs(NS(:,ii)) .* [rev==0]';
    NS_irrev( isnan( NS_irrev ) ) = 0; 
    if abs(sum(NS_irrev))~=sum(abs(NS_irrev))
       invalid_internal_cycles = invalid_internal_cycles + 1;
        % flux in the nullspace?
        disp(['Internal cycle #' num2str(ii) ' is invalid.'])
        C_pos = lab_flx(NS( :,ii )>0);
        C_neg = lab_flx(NS( :,ii )<0);
        disp(strjoin( C_pos , ' + '));
        disp([' - ' strjoin( C_neg , ' - ')]);
    else
        % flux in the nullspace?
        disp(['Internal cycle #' num2str(ii) ' is valid.'])
        C_pos = lab_flx(NS( :,ii )>0);
        C_neg = lab_flx(NS( :,ii )<0);
        disp(strjoin( C_pos , ' + '));
        disp([' - ' strjoin( C_neg , ' - ')]);
    end
end

disp(['Invalid internal cycles found: ' num2str(invalid_internal_cycles)])

if isempty(NS)
    disp('There are no internal cycles in this metabolic network')
end

end