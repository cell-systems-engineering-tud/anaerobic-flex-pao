function print_reactions_rev( N, v_lab, m_lab, rev )

% go through all columns and print the reactions

fprintf( 'Reactions of the network:\n');

for i=1:size( N, 2 ),   % i=1 ... first dimension of N = nr. of rows
    
    % positiv stoichiometric coeficent == substrate of the reaction
    p = find( N(:,i) > 0 );
    % positiv stoichiometric coeficent == product of the reaction
    s = find( N(:,i) < 0 );
    
    % we are at column i, print out the name of the reaction
    fprintf( '%.1d %s:\t', i, v_lab{i} );
    
    if ~isnan( s ),     % is there at least one substrate?
        % left side = substrates
        fprintf( '%.1d %s ', -N( s(1), i ), m_lab{ s(1) } ); % print first substrate
        
        for ii=2:length( s ) % are there more?
            fprintf( '+ %.1d %s ', -N( s(ii), i ), m_lab{ s(ii) } );
        end
    end
    
    if rev(i) > 0,
        fprintf( ' <--> ' );
    else
        fprintf( ' ---> ' );
    end
    
    if ~isnan( p ),     % is there at least one substrate?
        
        % right side = product
        fprintf( '%.1d %s ', N( p(1), i ), m_lab{ p(1) } ); % print first substrate
        
        for ii=2:length( p ) % are there more?
            fprintf( '+ %.1d %s ', N( p(ii), i ), m_lab{ p(ii) } );
        end
    end
    
    fprintf( '\n' );    % new line
    
end

if 0,        % put to 1 if flux and metabolite vector should be printed
    fprintf( '\nFlux vector\n' );    % new line
    
    for i=1:length(v_lab),
        fprintf( 1, '%.1d: %s\n', i, v_lab{i} )
    end
    
    fprintf( '\nMetabolite vector\n' );    % new line
    
    for i=1:length(m_lab),
        fprintf( 1, '%.1d: %s\n', i, m_lab{i} )
    end
end


end


