
runs = 15052:15097;
for i = 1:numel(runs)
    spe_file{i} = ['/home/dl11170/edatc/data/map' num2str(runs(i)) '_ei400.nxspe'];
end
par_file = '';
sqw_file = 'iron.sqw';
efix = 401;  % meV
emode = 1;   % 1=direct, 2=indirect
alatt = [2.87, 2.87, 2.87]; % in Angstrom
angdeg = [90, 90, 90];
u = [1,0,0];
v = [0,1,0];
psi = 0:2:90;
omega = 0; dpsi = 0; gl = 0; gs = 0;

gen_sqw (spe_file, par_file, sqw_file, efix, emode, alatt, angdeg,...
         u, v, psi, omega, dpsi, gl, gs)
     
%%
data_source = 'iron.sqw';

%proj = projaxes([1,0,0], [0,1,0], 'type', 'rrr');
%w1 = cut (data_source, proj, p1_bin, p2_bin, p3_bin, p4_bin)
%w2 = cut (data_source, p1_bin, p2_bin, p3_bin, p4_bin)
% w1, w2 equivalent because proj u and v same as u and v used in gen_sqw

% Here we use a new projection - plot axes along [110] and [-110]
% instead of [100] and [010].
proj = projaxes([1,1,0], [-1,1,0], 'type', 'rrr');
p1_bin = [-1.5,0.05,1.5];
p2_bin = [-1.5, 0.05, 1.5];
p3_bin = [-0.1, 0.1];
e_bin = [5,5,360];
w1 = cut_sqw(data_source, proj, p1_bin, p2_bin, p3_bin, e_bin)

plot(w1)

%%
w2 = cut(w1, [-2 0.05 2], [-2 0.05 2], [10 20])

%%
w3 = cut(w2, [], [100,120])
w4 = cut(w2, [], [120,140])
%%
acolor('black')
plot(w3)
%pl(w3)
acolor('r')
pd(w4)
legend({'w3','w4'})
title('')  % Remove default title
title('My Cool Figure!') % Add custom title
    
%%

w2 = cut(w1, [], [-1.1, -0.9], [])
plot(w2)

%%

proj = projaxes([1,1,0], [-1,1,0], 'type', 'rrr');
p1_bin = [-1.5,0.05,1.5];
p2_bin = [-1.5, 0.05, 1.5];
p3_bin = [-0.1, 0.1];
e_bin = [5,5,360];
w1dnd = cut_sqw(data_source, proj, p1_bin, p2_bin, p3_bin, e_bin, '-nopix')

plot(w1)

w2dnd = cut(w1dnd, [-2 0.025 2], [-2 0.05 2], [10 20])  % Ok because bin size the same

w2error = cut(w1dnd, [-2 0.025 2], [-2 0.05 2], [10 20])
% Error because w1dnd is dnd and we're changing the bin size of u1 axis.
