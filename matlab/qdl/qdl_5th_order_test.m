clear variables

dqmin = 0.01;
dqmax = 0.01;
maxerr = 0.001;

npt = 1e6;

%sys.add_node(i, C, G, H, p, B, q, S, dqmin, dqmax, maxerr);
%sys.add_branch(k, i, j, L, R, E, p, T, q, Z, dqmin, dqmax, maxerr);

%                                   Ra   La
%             .--------------------VVV---UUU-----.
%             |                       --->       |                wr               
%        .----o-------.                ia       +|        .-----.----.-----.
%        |    |1/Ro   |                     Ke*w( )       |     |    |     |
% Va/Ro (^)  [ ]     ===                     |  Kt*ia(^) Tm(v)  [ ]Bm === Jm
%        |    |   Clim|                          |        |     |    |     |
%        '----+-------'                                   '-----'----+-----'
%            ---                                                    ---
%             '                                                      '
%            Lb
%      o-----UUU----o      o 
%                    \     |
%                         ===   --->   
%                          |
%                         ---
%                          '

Va = 10.0;
Ro = 0.001;
Clim = 0.01;

Ra = 0.5;
La = 0.1;

Tm = 1.0;
Bm = 0.01;
Jm = 0.1;
Kt = 0.1;
Ke = 0.1;


sys = QdlSystem(2, 1, npt);

                   % k  i  j  a     b     c   dk  d   ek  e   src min  max  err
idx = sys.add_node  (1,       0,    0,    0,  0,  0,  0,  0,  10, 0.1, 0.1, maxerr);
idx = sys.add_branch(1, 1, 0, La,   Ra,   0,  2, -Ke, 0,  0,   0, 0.1, 0.1, maxerr);
idx = sys.add_node  (2,       Jm,   Bm,   0,  1,  0,  1,  Kt,  0, 0.1, 0.1, maxerr);

sys.swatom(1) = 1;
sys.swper(1)  = 1e-2;
sys.swduty(1) = 0.5;
sys.swon(1)   = 10.0;
sys.swoff(1)  = 1.0;

sys.init();
sys.runto(0.5);
%sys.H(2) = -1;
%sys.runto(100);
%sys.source(1) = 20.0;
%sys.runto(150);

for k=1:sys.n
    plot(sys.tout(k,1:sys.iout(k)), sys.qout(k,1:sys.iout(k)), 'b.-'); hold on;
end

figure;
for k=1:sys.n
    plot(sys.tupd(k, 1:sys.iupd(k)), cumsum(sys.nupd(k, 1:sys.iupd(k)))); hold on;
end
