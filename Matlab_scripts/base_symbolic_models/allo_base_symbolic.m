% for allos with 2 HEs, classification of fixed points using linear stability
% analysis, the Jacobian matrix, and eigenvectors

% this code has not been heavily optimized and may be slow
% prints solutions and their stabilities for each s value

s_val_range = 1e-3; % scalar or vector (e.g. logspace) of s values

mu_val = 1e-8; % constant value of forward mutation rate
nu_val = 1e-8; % constant value of backward mutation rate

h1_val = 0; % h1 dominance coefficient value, constant
h2_val = 0; % h2 dominance coefficient value, constant
h3_val = 0; % h3 dominance coefficient value, constant

syms g00 g01 g10 g11 s h1 h2 h3 mu nu

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
sel_g00 = g00^2*w0+g00*g01*w1+g00*g10*w1+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2;
sel_g10 = g00*g10*w1+g10^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g10*g11*w3;
sel_g01 = g00*g01*w1+g01^2*w2+(1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g01*g11*w3;
sel_g11 = (1/2)*g01*g10*w2+(1/2)*g00*g11*w2+g01*g11*w3+g10*g11*w3+g11^2*w4;

% equations for mutation
mut_g00 = sel_g00*(1-mu)^2 + sel_g01*(1-mu)*nu + sel_g10*(1-mu)*nu + sel_g11*nu^2 - g00;
mut_g01 = sel_g00*mu*(1-mu) + sel_g01*(1-mu)*(1-nu) + sel_g10*mu*nu + sel_g11*(1-nu)*nu - g01;
mut_g10 = sel_g00*mu*(1-mu) + sel_g01*mu*nu + sel_g10*(1-mu)*(1-nu) + sel_g11*(1-nu)*nu - g10;
mut_g11 = sel_g00*mu^2 + sel_g01*mu*(1-nu) + sel_g10*mu*(1-nu) + sel_g11*(1-nu)^2 - g11;

mut_eqn_set = [mut_g00, mut_g01, mut_g10, mut_g11];

for i = 1:length(mut_eqn_set)
    % removes g10 from the equation by replacing it with 1-(g00+g01+g11)
    mut_eqn_set(i) = subs(mut_eqn_set(i), g10, 1-(g00+g01+g11));    
end


%creates the Jacobian of the system
jac_matrix = [diff(mut_eqn_set(1), g00), diff(mut_eqn_set(1), g01), diff(mut_eqn_set(1), g11); 
                diff(mut_eqn_set(2), g00), diff(mut_eqn_set(2), g01), diff(mut_eqn_set(2), g11); 
                diff(mut_eqn_set(4), g00), diff(mut_eqn_set(4), g01), diff(mut_eqn_set(4), g11)];

neutral_stable_g00 = [];
neutral_stable_g01 = [];
neutral_stable_g11 = [];
neutral_stable_s = [];

selected_stable_g00 = [];
selected_stable_g01 = [];
selected_stable_g11 = [];
selected_stable_s = [];

unstable_g00 = [];
unstable_g01 = [];
unstable_g11 = [];
unstable_s = [];

for i = 1:length(s_val_range)

    %solves for the fixed points of the system
    [g00_root_vals, g01_root_vals, g11_root_vals] = root_solns(mut_eqn_set(1), mut_eqn_set(2), mut_eqn_set(4), mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, g00, g01, g11)
        
    %evaluating the jacobian and stability of each fixed point
    [fixed_pt_stabilities] = linear_stability_analysis(jac_matrix, mu, mu_val, nu, nu_val, s, s_val_range(i), h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals) 

    for j = 1:length(fixed_pt_stabilities)

        if fixed_pt_stabilities(j) == 0
            unstable_g00(end+1) = g00_root_vals(j);
            unstable_g01(end+1) = g01_root_vals(j);
            unstable_g11(end+1) = g11_root_vals(j);
            unstable_s(end+1) = s_val_range(i);

        elseif fixed_pt_stabilities(j) == 1
            if g00_root_vals(j) > .3333
                selected_stable_g00(end+1) = g00_root_vals(j);
                selected_stable_g01(end+1) = g01_root_vals(j);
                selected_stable_g11(end+1) = g11_root_vals(j);
                selected_stable_s(end+1) = s_val_range(i);
            else
                neutral_stable_g00(end+1) = g00_root_vals(j);
                neutral_stable_g01(end+1) = g01_root_vals(j);
                neutral_stable_g11(end+1) = g11_root_vals(j);
                neutral_stable_s(end+1) = s_val_range(i);
            end
        end
    end
end

% figure
% 
% plot(neutral_stable_s, neutral_stable_g00+neutral_stable_g01)
% hold on
% plot(selected_stable_s, selected_stable_g00+selected_stable_g01)
% plot(unstable_s, unstable_g00+unstable_g01, 'LineStyle','--')
% 
% xscale log
% 
% figure
% 
% plot(neutral_stable_s, neutral_stable_g00)
% hold on
% plot(selected_stable_s, selected_stable_g00)
% plot(unstable_s, unstable_g00, 'LineStyle','--')
% 
% xscale log
% 
% figure
% 
% plot(neutral_stable_s, neutral_stable_g01)
% hold on
% plot(selected_stable_s, selected_stable_g01)
% plot(unstable_s, unstable_g01, 'LineStyle','--')
% 
% xscale log
% 
% figure
% 
% plot(neutral_stable_s, neutral_stable_g11)
% hold on
% plot(selected_stable_s, selected_stable_g11)
% plot(unstable_s, unstable_g11, 'LineStyle','--')
% 
% xscale log



%%% FUNCTIONS %%%

function [g00_root_vals, g01_root_vals, g11_root_vals] = root_solns(mut_g00_eqn, mut_g01_eqn, mut_g11_eqn, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, g00, g01, g11)

    %function which uses vpasolve to find the fixed points/roots of the system

    g00_eqn = subs(mut_g00_eqn, mu, mu_val);
    g00_eqn = subs(g00_eqn, nu, nu_val);
    g00_eqn = subs(g00_eqn, s, s_val);
    g00_eqn = subs(g00_eqn, h1, h1_val);
    g00_eqn = subs(g00_eqn, h2, h2_val);
    g00_eqn = subs(g00_eqn, h3, h3_val);

    g01_eqn = subs(mut_g01_eqn, mu, mu_val);
    g01_eqn = subs(g01_eqn, nu, nu_val);
    g01_eqn = subs(g01_eqn, s, s_val);
    g01_eqn = subs(g01_eqn, h1, h1_val);
    g01_eqn = subs(g01_eqn, h2, h2_val);
    g01_eqn = subs(g01_eqn, h3, h3_val);

    g11_eqn = subs(mut_g11_eqn, mu, mu_val);
    g11_eqn = subs(g11_eqn, nu, nu_val);
    g11_eqn = subs(g11_eqn, s, s_val);
    g11_eqn = subs(g11_eqn, h1, h1_val);
    g11_eqn = subs(g11_eqn, h2, h2_val);
    g11_eqn = subs(g11_eqn, h3, h3_val);


    [g00_root_vals, g01_root_vals, g11_root_vals] = vpasolve([g00_eqn, g01_eqn, g11_eqn], [g00, g01, g11]);
end


function [pd_value] = pd_evaluation(jacobian_entry, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_root_val, g01, g01_root_val, g11, g11_root_val)
    
    %%%function which evaluates a partial derivative by substituting a root of
    %the system
    %used to evaluate the jacobian at one entry%%%

    pd_value = subs(jacobian_entry, mu, mu_val);
    pd_value = subs(pd_value, nu, nu_val);
    pd_value = subs(pd_value, s, s_val);
    pd_value = subs(pd_value, h1, h1_val);
    pd_value = subs(pd_value, h2, h2_val);
    pd_value = subs(pd_value, h3, h3_val);
    pd_value = subs(pd_value, g00, g00_root_val);
    pd_value = subs(pd_value, g01, g01_root_val);
    pd_value = subs(pd_value, g11, g11_root_val);
end

function [fixed_pt_stabilities] = linear_stability_analysis(jacobian_matrix, mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_root_vals, g01, g01_root_vals, g11, g11_root_vals)
    
    %%%evaluates the jacobian in full by calling pd_evaluation
    %Then, calculates the eigenvalues and vectors of the Jacobian
    %Using the eigenvalues, classifies stability
    %Returns a 0 if the fixed point is unstable, a 1 if stable

    fixed_pt_stabilities = g00_root_vals;

    jacobian_eval = [0, 0, 0; 
                     0, 0, 0
                     0, 0, 0];

    for i = 1:length(g00_root_vals)
        for j = 1:length(jacobian_eval)
            for k = 1:length(jacobian_eval)
                jacobian_eval(j, k) = pd_evaluation(jacobian_matrix(j, k), mu, mu_val, nu, nu_val, s, s_val, h1, h1_val, h2, h2_val, h3, h3_val, g00, g00_root_vals(i), g01, g01_root_vals(i), g11, g11_root_vals(i)); 
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
