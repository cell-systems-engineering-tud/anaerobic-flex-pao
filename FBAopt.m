function [v_opt, Aeq] = FBAopt(N,objective,meas_idx,meas_flx,OFF_idx,rev,upperbounded_idx,upperbounded_flx)

% perform optimization
f = zeros( size( N, 2 ), 1 );  %creates a vertical vector the same number of lines as matrix N columns zeros
f(abs(objective)) = 1 * abs(objective)/objective;    %Linear objective function vector f, mu needs to be maximized

Ceq = eye(size(N,2)); %same size as the number of fluxes
M = Ceq(meas_idx,:); %measurements
turned_off_fluxes = Ceq(OFF_idx,:); %knocked-out fluxes are turned off

Ceq = [M;turned_off_fluxes];

Aeq = [ N; Ceq ]; %Adds the new line to N      %Matrix for linear equality constraints
beq = [ zeros( size( N,1 ), 1 ); meas_flx ; zeros(length(OFF_idx),1)]; % respect stoichiometry and set maintenance flx  %Vector for linear equality constraints

n_rev = find(rev); %says in which columns of rev there are ones (reversibilities)
n_irr = setdiff(1:length(rev),n_rev);
Nirr = eye( size(N, 2) );    %Matrix for linear inequality constraints
Nirr = -Nirr( n_irr, : );

Cie = eye( size(N,2) ); %same size as the number of fluxes
Cie_upper = Cie( upperbounded_idx, : ); %inequality constraints
Aie = [ Nirr; Cie_upper];
bie = [ zeros( length(n_irr), 1 ); upperbounded_flx]; % Aie*x <= bie;  %Vector for linear inequality constraints

v_opt = linprog( -f, Aie, bie, Aeq, beq);
end