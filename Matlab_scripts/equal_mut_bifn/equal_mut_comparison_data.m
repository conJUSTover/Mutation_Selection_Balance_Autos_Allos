s_val_range = logspace(-9, -3, 500);
mu_val = 1e-8;
nu_val = 1e-8;

% recessive case
h_val = 0;
h1_val = 0;
h2_val = 0;
h3_val = 0;

dip_rec = diploids_bifn_data_stable_only(s_val_range, mu_val, nu_val, h_val);
disp("Diploid recessive complete")
auto_rec = auto_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp("Auto recessive complete")
allo_rec = allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp("Allo recessive complete")

writematrix(dip_rec, 'dip_rec.csv')
writematrix(auto_rec, 'auto_rec.csv')
writematrix(allo_rec, 'allo_rec.csv')

disp('Recessive case complete')

% additive case
h_val = .5;
h1_val = .25;
h2_val = .5;
h3_val = .75;

dip_add = diploids_bifn_data_stable_only(s_val_range, mu_val, nu_val, h_val);
auto_add = auto_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
allo_add = allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

writematrix(dip_add, 'dip_add.csv')
writematrix(auto_add, 'auto_add.csv')
writematrix(allo_add, 'allo_add.csv')

disp('Additive case complete')

% dominant case
h_val = 1;
h1_val = 1;
h2_val = 1;
h3_val = 1;

dip_dom = diploids_bifn_data_stable_only(s_val_range, mu_val, nu_val, h_val);
auto_dom = auto_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
allo_dom = allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);

writematrix(dip_dom, 'dip_dom.csv')
writematrix(auto_dom, 'auto_dom.csv')
writematrix(allo_dom, 'allo_dom.csv')

