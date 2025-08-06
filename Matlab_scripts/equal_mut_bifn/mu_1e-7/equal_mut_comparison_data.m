s_val_range = logspace(-9, -3, 500);
mu_val = 1e-7;
nu_val = 1e-7;

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

% partially recessive case
h_val = .4;
h1_val = .2;
h2_val = .4;
h3_val = .6;

dip_part_rec = diploids_bifn_data_stable_only(s_val_range, mu_val, nu_val, h_val);
disp("Diploid partially recessive complete")
auto_part_rec = auto_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp("Auto partially recessive complete")
allo_part_rec = allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp("Allo partially recessive complete")

writematrix(dip_part_rec, 'dip_part_rec.csv')
writematrix(auto_part_rec, 'auto_part_rec.csv')
writematrix(allo_part_rec, 'allo_part_rec.csv')

% additive case
h_val = .5;
h1_val = .25;
h2_val = .5;
h3_val = .75;

dip_add = diploids_bifn_data_stable_only(s_val_range, mu_val, nu_val, h_val);
disp("Diploid additive complete")
auto_add = auto_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp("Auto additive complete")
allo_add = allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp("Allo additive complete")

writematrix(dip_add, 'dip_add.csv')
writematrix(auto_add, 'auto_add.csv')
writematrix(allo_add, 'allo_add.csv')

% partially dominant case
h_val = .6;
h1_val = .4;
h2_val = .6;
h3_val = .8;

dip_part_dom = diploids_bifn_data_stable_only(s_val_range, mu_val, nu_val, h_val);
disp("Diploid partially dominant complete")
auto_part_dom = auto_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp("Auto partially dominant complete")
allo_part_dom = allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp("Allo partially dominant complete")

writematrix(dip_part_dom, 'dip_part_dom.csv')
writematrix(auto_part_dom, 'auto_part_dom.csv')
writematrix(allo_part_dom, 'allo_part_dom.csv')

% dominant case
h_val = 1;
h1_val = 1;
h2_val = 1;
h3_val = 1;

dip_dom = diploids_bifn_data_stable_only(s_val_range, mu_val, nu_val, h_val);
disp("Diploid dominant complete")
auto_dom = auto_bifn_data_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp("Auto dominant complete")
allo_dom = allo_bifn_stable_only(s_val_range, mu_val, nu_val, h1_val, h2_val, h3_val);
disp("Allo dominant complete")

writematrix(dip_dom, 'dip_dom.csv')
writematrix(auto_dom, 'auto_dom.csv')
writematrix(allo_dom, 'allo_dom.csv')

