function [neutral_data, selected_data, unstable_data] = auto_bifn_data(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val)
%Generates autopolyploid bifurcation plotting data
%   for autos, classification of fixed points using linear stability
%   analysis, the Jacobian matrix, and eigenvectors

g0 = sym('g0');
g1 = sym('g1');
g2 = sym('g2');
G0 = sym('G0');
G1 = sym('G1');
G2 = sym('G2');
G3 = sym('G3');
G4 = sym('G4');
s = sym('s');
h1 = sym('h1');
h2 = sym('h2');
h3 = sym('h3');
mu = sym('mu');
nu = sym('nu');

% assumptions on the parameters of the model; theoretical bounds
assume(g0>=0 & g0<=1);
assume(g1>=0 & g1<=1);
assume(g2>=0 & g2<=1);
assume(s>=-1 & s<=1);
assume(h1>=0 & h1<=1);
assume(h2>=0 & h2<=1);
assume(h3>=0 & h3<=1);
assume(mu>=0 & mu<=1);
assume(nu>=0 & nu<=1);
assume(G0>=0 & G0<=1);
assume(G1>=0 & G1<=1);
assume(G2>=0 & G2<=1);
assume(G3>=0 & G3<=1);
assume(G4>=0 & G4<=1);

% equations to parameterize relative fitnesses
wbar = 1 - s*(G1*h1 + G2*h2 + G3*h3 + G4);
w0 = 1/wbar;
w1 = (1-s*h1)/wbar;
w2 = (1-s*h2)/wbar;
w3 = (1-s*h3)/wbar;
w4 = (1-s)/wbar;

% equations for selection
sel_g0 = G0*w0 + (1/2)*G1*w1 + (1/6)*G2*w2;
sel_g1 = (1/2)*G1*w1 + (2/3)*G2*w2 + (1/2)*G3*w3;
sel_g2 = (1/6)*G2*w2 + (1/2)*G3*w3 + G4*w4;

% equations for mutation
mut_g0 = sel_g0*((1-mu)^2) + sel_g1*(1-mu)*nu + sel_g2*(nu^2) - g0;
mut_g1 = 2*sel_g0*(1-mu)*mu + sel_g1*(1-mu)*(1-nu)+2*sel_g2*(1-nu)*nu - g1;
mut_g2 = sel_g0*(mu^2) + sel_g1*mu*(1-nu) + sel_g2*((1-nu)^2) - g2;

mut_exp_set = [mut_g0, mut_g1, mut_g2];

%substituing genotypes for gametes and removing g2 using g0+g1+g2 = 1
for i = 1:length(mut_exp_set)
    mut_exp_set(i) = subs(mut_exp_set(i), G0, g0^2);
    mut_exp_set(i) = subs(mut_exp_set(i), G1, 2*g0*g1);
    mut_exp_set(i) = subs(mut_exp_set(i), G2, (2*g0*g2 + g1^2));
    mut_exp_set(i) = subs(mut_exp_set(i), G3, 2*g1*g2);
    mut_exp_set(i) = subs(mut_exp_set(i), G4, g2^2);
    mut_exp_set(i) = subs(mut_exp_set(i), g2, (1-g1-g0));
end


%creates the Jacobian of the system
jac_matrix = [diff(mut_exp_set(1), g0), diff(mut_exp_set(1), g1); diff(mut_exp_set(2), g0), diff(mut_exp_set(2), g1)];

%creates a partially evaluated vector of ODEs for constant parameters
mut_exp_set_eval = mut_exp_set;

for i = 1:length(mut_exp_set_eval)
    mut_exp_set_eval(i) = auto_diff_eq_eval(mut_exp_set(i), mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val);
end

%similarly, creates a partially evaluted Jacobian for constant parameters
jac_matrix_eval = jac_matrix;

for i = 1:length(jac_matrix_eval)
    for j = 1:length(jac_matrix_eval)
        jac_matrix_eval(i, j) = auto_pd_evaluation_init(jac_matrix_eval(i,j), mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val);
    end
end

neutral_g0 = [];
neutral_g1 = [];
neutral_s = [];

selected_g0 = [];
selected_g1 = [];
selected_s = [];

unstable_g0 = [];
unstable_g1 = [];
unstable_s = [];

for i = 1:length(s_val_range)

    %solves for the fixed points of the system
    [g0_root_vals, g1_root_vals] = auto_root_solns(mut_exp_set_eval(1), mut_exp_set_eval(2), s, s_val_range(i), g0, g1);
        
    %evaluating the jacobian and stability of each fixed point
    [fixed_pt_stabilities] = auto_linear_stability_analysis(jac_matrix_eval, s, s_val_range(i), g0, g0_root_vals, g1, g1_root_vals); 

    for j = 1:length(fixed_pt_stabilities)

        if fixed_pt_stabilities(j) == 0
            unstable_g0(end+1) = g0_root_vals(j);
            unstable_g1(end+1) = g1_root_vals(j);
            unstable_s(end+1) = s_val_range(i);

        elseif fixed_pt_stabilities(j) == 1
            if g0_root_vals(j) > 1/3 && s_val_range(i) > 0
                selected_g0(end+1) = g0_root_vals(j);
                selected_g1(end+1) = g1_root_vals(j);
                selected_s(end+1) = s_val_range(i);
            elseif g0_root_vals(j) > 2/3 && s_val_range(i) < 0
                selected_g0(end+1) = g0_root_vals(j);
                selected_g1(end+1) = g1_root_vals(j);
                selected_s(end+1) = s_val_range(i);
            else
                neutral_g0(end+1) = g0_root_vals(j);
                neutral_g1(end+1) = g1_root_vals(j);
                neutral_s(end+1) = s_val_range(i);
            end
        end
    end
end

neutral_q = ones(1, length(neutral_g0)) - neutral_g0 - .5*neutral_g1;
selected_q = ones(1, length(selected_g0)) - selected_g0 - .5*selected_g1;
unstable_q = ones(1, length(unstable_g0)) - unstable_g0 - .5*unstable_g1;

neutral_g2 = 1 - neutral_g0 - neutral_g1;
selected_g2 = 1 - selected_g0 - selected_g1;
unstable_g2 = 1 - unstable_g0 - unstable_g1;

neutral_genotypes = [neutral_g0.^2; 2*neutral_g0.*neutral_g1; neutral_g1.^2+2*neutral_g0.*neutral_g2; 2*neutral_g1.*neutral_g2; neutral_g2.^2];
selected_genotypes = [selected_g0.^2; 2*selected_g0.*selected_g1; selected_g1.^2+2*selected_g0.*selected_g2; 2*selected_g1.*selected_g2; selected_g2.^2];
unstable_genotypes = [unstable_g0.^2; 2*unstable_g0.*unstable_g1; unstable_g1.^2+2*unstable_g0.*unstable_g2; 2*unstable_g1.*unstable_g2; unstable_g2.^2];

neutral_abs_fitness = [ones(1, length(neutral_s)); 1-h1_val*neutral_s; 1-h2_val*neutral_s; 1-h3_val*neutral_s; 1-neutral_s];
selected_abs_fitness = [ones(1, length(selected_s)); 1-h1_val*selected_s; 1-h2_val*selected_s; 1-h3_val*selected_s; 1-selected_s];
unstable_abs_fitness = [ones(1, length(unstable_s)); 1-h1_val*unstable_s; 1-h2_val*unstable_s; 1-h3_val*unstable_s; 1-unstable_s];

neutral_w_bar = zeros(1, length(neutral_s));
selected_w_bar = zeros(1, length(selected_s));
unstable_w_bar = zeros(1, length(unstable_s));

for i = 1:length(neutral_w_bar)
    neutral_w_bar(i) = dot(neutral_genotypes(:, i), neutral_abs_fitness(:, i));
end

for i = 1:length(selected_w_bar)
    selected_w_bar(i) = dot(selected_genotypes(:, i), selected_abs_fitness(:, i));
end

for i = 1:length(unstable_w_bar)
    unstable_w_bar(i) = dot(unstable_genotypes(:, i), unstable_abs_fitness(:, i));
end

neutral_data_transposed = vertcat(neutral_s, neutral_q, neutral_g0, neutral_g1, neutral_g2, neutral_w_bar);
neutral_data = neutral_data_transposed.';

selected_data_transposed = vertcat(selected_s, selected_q, selected_g0, selected_g1, selected_g2, selected_w_bar);
selected_data = selected_data_transposed.';

unstable_data_transposed = vertcat(unstable_s, unstable_q, unstable_g0, unstable_g1, unstable_g2, unstable_w_bar);
unstable_data = unstable_data_transposed.';

end



%%% FUNCTIONS %%%

function [diff_eqn_eval] = auto_diff_eq_eval(diff_eqn, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val)

    % substitutes all constant parameters to create simplified differential
    % equations

    diff_eqn_eval = subs(diff_eqn, mu, mu_val);
    diff_eqn_eval = subs(diff_eqn_eval, nu, nu_val);
    diff_eqn_eval = subs(diff_eqn_eval, h1, h1_val);
    diff_eqn_eval = subs(diff_eqn_eval, h2, h2_val);
    diff_eqn_eval = subs(diff_eqn_eval, h3, h3_val);

end

function [g0_root_vals, g1_root_vals] = auto_root_solns(g0_eqn_eval, g1_eqn_eval, s, s_val, g0, g1)

    %function which uses vpasolve to find the fixed points/roots of the system

    g0_eqn = subs(g0_eqn_eval, s, s_val);

    g1_eqn = subs(g1_eqn_eval, s, s_val);

    [g0_root_vals, g1_root_vals] = vpasolve([g0_eqn, g1_eqn], [g0, g1]);
end


function [pd_value] = auto_pd_evaluation_init(jacobian_entry, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, mu, mu_val);
    pd_value = subs(pd_value, nu, nu_val);
    pd_value = subs(pd_value, h1, h1_val);
    pd_value = subs(pd_value, h2, h2_val);
    pd_value = subs(pd_value, h3, h3_val);

end

function [pd_value] = auto_pd_evaluation_final(jacobian_entry, s, s_val, g0, g0_root_val, g1, g1_root_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, s, s_val);
    pd_value = subs(pd_value, g0, g0_root_val);
    pd_value = subs(pd_value, g1, g1_root_val);
end

function [fixed_pt_stabilities] = auto_linear_stability_analysis(jacobian_matrix, s, s_val, g0, g0_root_vals, g1, g1_root_vals)
    
    %%%evaluates the jacobian in full by calling pd_evaluation
    %Then, calculates the eigenvalues and vectors of the Jacobian
    %Using the eigenvalues, det(J) and tr(J), classifies stability
    %Returns a 0 if the fixed point is unstable, a 1 if stable

    fixed_pt_stabilities = g0_root_vals;

    jacobian_eval = [0, 0; 0, 0];

    for i = 1:length(g0_root_vals)
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
                jacobian_eval(j, k) = auto_pd_evaluation_final(jacobian_matrix(j, k), s, s_val, g0, g0_root_vals(i), g1, g1_root_vals(i)); 
            end
        end

        %calulating the trace and determinant of the evaluated jacobian
        trace_jac = trace(jacobian_eval);
        det_jac = det(jacobian_eval);

        %classifies the fixed point according to the trace and determinant
        if det_jac < 0
            fixed_pt_stabilities(i) = 0; %0 for unstable saddle point
        elseif det_jac > 0
            if trace_jac < 0 && trace_jac^2 - 4*det_jac > 0
                fixed_pt_stabilities(i) = 1; %1 for stable node
            else
                disp('Error. Unexpected stability type from linear stability analysis.')
            end
        else
            disp('Error. Unexpected stability type from linear stability analysis.')
        end
    end
end
