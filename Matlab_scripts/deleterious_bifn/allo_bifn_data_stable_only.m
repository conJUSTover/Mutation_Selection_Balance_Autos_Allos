function stable_data = allo_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val)
%Generates allo bifurcation plotting data
%   classifies fixed points using linear stability analysis, the Jacobian 
%   matrix, and eigenvectors

g00 = sym('g00'); % gamete frequency with 0 derived alleles
g01 = sym('g01'); % gamete frequency with 1 derived allele in 'b' subgenome
g10 = sym('g10'); % gamete frequency with 1 derived allele in 'a' subgenome
g11 = sym('g11'); % gamete frequency with 2 derived alleles
s = sym('s'); % selection coefficient
h1 = sym('h1'); % dominance coefficient for genotype with 1 derived allele
h2 = sym('h2'); % dominance coefficient for genotype with 2 derived alleles
h3 = sym('h3'); % dominance coefficient for genotype with 3 derived alleles
mu = sym('mu'); % forward mutation rate
nu = sym('nu'); % backward mutation rate

% assumptions on the parameters of the model; theoretical bounds
assume(g00>=0 & g00<=1);
assume(g10>=0 & g10<=1);
assume(g01>=0 & g01<=1);
assume(g11>=0 & g11<=1);
assume(s>=0 & s<=1);
assume(h1>=0 & h1<=1);
assume(h2>=0 & h2<=1);
assume(h3>=0 & h3<=1);
assume(mu>=0 & mu<=1);

% equations to parameterize relative fitnesses
wbar = (1-2*s*(h1*(g00*g10+g00*g01)+h2*(g00*g11+g01*g10)+h3*(g01*g11+g10*g11))-s*(h2*(g01^2+g10^2)+g11^2));
w0 = 1/wbar;
w1 = (1-s*h1)/wbar;
w2 = (1-s*h2)/wbar;
w3 = (1-s*h3)/wbar;
w4 = (1-s)/wbar;

% equations for selection
sel_g00 = g00^2*w0 + (.5)*2*g01*g00*w1 + (.5)*2*g10*g00*w1 + (.25)*(2*g00*g11+2*g01*g10)*w2;
sel_g10 = (.5)*2*g10*g00*w1 + g10^2*w2 + (.25)*(2*g00*g11+2*g01*g10)*w2 + (.5)*2*g10*g11*w3;
sel_g01 = (.5)*2*g01*g00*w1 + g01^2*w2 + (.25)*(2*g00*g11+2*g01*g10)*w2 + (.5)*2*g01*g11*w3;
sel_g11 = (.25)*(2*g00*g11+2*g01*g10)*w2 + (.5)*2*g01*g11*w3 + (.5)*2*g10*g11*w3 + g11^2*w4;

% equations for mutation
mut_g00 = sel_g00*(1-mu)^2 + sel_g01*(1-mu)*nu + sel_g10*(1-mu)*nu + sel_g11*nu^2 - g00;
mut_g01 = sel_g00*mu*(1-mu) + sel_g01*(1-mu)*(1-nu) + sel_g10*mu*nu + sel_g11*(1-nu)*nu - g01;
mut_g10 = sel_g00*mu*(1-mu) + sel_g01*mu*nu + sel_g10*(1-mu)*(1-nu) + sel_g11*(1-nu)*nu - g10;
mut_g11 = sel_g00*mu^2 + sel_g01*mu*(1-nu) + sel_g10*mu*(1-nu) + sel_g11*(1-nu)^2 - g11;

mut_exp_set = [mut_g00, mut_g01, mut_g10, mut_g11];

for i = 1:length(mut_exp_set)
    % removes g10 from the equation by replacing it with 1-(g00+g01+g11)
    mut_exp_set(i) = subs(mut_exp_set(i), g10, 1-(g00+g01+g11));    
end


%creates the Jacobian of the system
jac_matrix = [diff(mut_exp_set(1), g00), diff(mut_exp_set(1), g01), diff(mut_exp_set(1), g11); 
                diff(mut_exp_set(2), g00), diff(mut_exp_set(2), g01), diff(mut_exp_set(2), g11); 
                diff(mut_exp_set(4), g00), diff(mut_exp_set(4), g01), diff(mut_exp_set(4), g11)];

%creates a partially evaluated vector of ODEs for constant parameters
mut_exp_set_eval = mut_exp_set;

for i = 1:length(mut_exp_set_eval)
    mut_exp_set_eval(i) = allo_diff_eq_eval(mut_exp_set(i), mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val);
end

%similarly, creates a partially evaluted Jacobian for constant parameters
jac_matrix_eval = jac_matrix;

for i = 1:length(jac_matrix_eval)
    for j = 1:length(jac_matrix_eval)
        jac_matrix_eval(i, j) = allo_pd_evaluation_init(jac_matrix_eval(i,j), mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val);
    end
end

stable_g00 = zeros(1, length(s_val_range));
stable_g01 = zeros(1, length(s_val_range));
stable_g11 = zeros(1, length(s_val_range));
stable_s = s_val_range;

for i = 1:length(s_val_range)

    %solves for the fixed points of the system
    [g00_root_vals, g01_root_vals, g11_root_vals] = allo_root_solns(mut_exp_set_eval(1), mut_exp_set_eval(2), mut_exp_set_eval(4), s, s_val_range(i), g00, g01, g11);
        
    if length(g00_root_vals) > 1
        disp('Error. Unexpected number of fixed points found.');
        return;
    end

    stable_g00(1, i) = g00_root_vals;
    stable_g01(1, i) = g01_root_vals;
    stable_g11(1, i) = g11_root_vals;
end

stable_q = stable_g11 + stable_g01;

stable_g10 = 1 - stable_g00 - stable_g01 - stable_g11;

stable_G1 = 2*stable_g00.*stable_g10 + 2*stable_g00.*stable_g01;
stable_G2 = stable_g10.^2 + 2*(stable_g00.*stable_g11 + stable_g01.*stable_g10) + stable_g01.^2;
stable_G3 = 2*stable_g11.*stable_g10 + 2*stable_g11.*stable_g01;
stable_genotypes = [stable_g00.^2; stable_G1; stable_G2; stable_G3; stable_g11.^2];

stable_abs_fitness = [ones(1, length(stable_s)); 1-h1_val*stable_s; 1-h2_val*stable_s; 1-h3_val*stable_s; 1-stable_s];

stable_w_bar = zeros(1, length(stable_s));

for i = 1:length(stable_w_bar)
    stable_w_bar(i) = dot(stable_genotypes(:, i), stable_abs_fitness(:, i));
end

stable_data_transposed = vertcat(stable_s, stable_q, stable_g00, stable_g01, stable_g10, stable_g11, stable_w_bar);
stable_data = stable_data_transposed.';

end


%%% FUNCTIONS %%%

function [diff_eqn_eval] = allo_diff_eq_eval(diff_eqn, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val)

    % substitutes all constant parameters to create simplified differential
    % equations

    diff_eqn_eval = subs(diff_eqn, mu, mu_val);
    diff_eqn_eval = subs(diff_eqn_eval, nu, nu_val);
    diff_eqn_eval = subs(diff_eqn_eval, h1, h1_val);
    diff_eqn_eval = subs(diff_eqn_eval, h2, h2_val);
    diff_eqn_eval = subs(diff_eqn_eval, h3, h3_val);

end

function [g00_root_vals, g01_root_vals, g11_root_vals] = allo_root_solns(mut_g00_eqn, mut_g01_eqn, mut_g11_eqn, s, s_val, g00, g01, g11)

    %function which uses vpasolve to find the fixed points/roots of the system

    g00_eqn = subs(mut_g00_eqn, s, s_val);

    g01_eqn = subs(mut_g01_eqn, s, s_val);

    g11_eqn = subs(mut_g11_eqn, s, s_val);

    [g00_root_vals, g01_root_vals, g11_root_vals] = vpasolve([g00_eqn, g01_eqn, g11_eqn], [g00, g01, g11]);
end


function [pd_value] = allo_pd_evaluation_init(jacobian_entry, mu, mu_val, nu, nu_val, h1, h1_val, h2, h2_val, h3, h3_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, mu, mu_val);
    pd_value = subs(pd_value, nu, nu_val);
    pd_value = subs(pd_value, h1, h1_val);
    pd_value = subs(pd_value, h2, h2_val);
    pd_value = subs(pd_value, h3, h3_val);

end

function [pd_value] = allo_pd_evaluation_final(jacobian_entry, s, s_val, g00, g00_root_val, g01, g01_root_val, g11, g11_root_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, s, s_val);
    pd_value = subs(pd_value, g00, g00_root_val);
    pd_value = subs(pd_value, g01, g01_root_val);
    pd_value = subs(pd_value, g11, g11_root_val);
end

function [fixed_pt_stabilities] = allo_linear_stability_analysis(jacobian_matrix, s, s_val, g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals)
    
    %%%evaluates the jacobian in full by calling pd_evaluation
    %Then, calculates the eigenvalues and vectors of the Jacobian
    %Using the eigenvalues, classifies stability
    %Returns a 0 if the fixed point is unstable, a 1 if stable

    fixed_pt_stabilities = g00_root_vals;

    jacobian_eval = [0, 0, 0; 
                     0, 0, 0;
                     0, 0, 0];

    for i = 1:length(g00_root_vals)
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
                jacobian_eval(j, k) = allo_pd_evaluation_final(jacobian_matrix(j, k), s, s_val, g00, g00_root_vals(i), g01, g01_root_vals(i), g11, g11_root_vals(i)); 
            end
        end

        %calulating the eigenvalues of the evaluated jacobian
        [eigenvectors, eigenvalues] = eig(jacobian_eval);

        %classifies the fixed point according to the trace and determinant
        if eigenvalues(1,1) < 0 && eigenvalues(2,2) < 0 && eigenvalues(3,3) < 0
            fixed_pt_stabilities(i) = 1; %0 for stable node
        elseif eigenvalues(1,1) > 0 && eigenvalues(2,2) > 0 && eigenvalues(3,3) > 0
            disp('Error. Unexpected stability type from linear stability analysis.')
        else
            fixed_pt_stabilities(i) = 0; %0 for unstable saddle point
        end
    end
end
