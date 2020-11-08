% System params
k = 0.1;
tau = 0.05;
dt = 0.01;
ns = 6;
nu = 2;
N = 150;

% Matrices
A_d = eye(ns)+dt*[0 0 1 0 0 0;0 0 0 1 0 0;0 0 -k 0 1 0;0 0 0 -k 0 1;0 0 0 0 -1/tau 0; 0 0 0 0 0 -1/tau];
B_d = zeros(ns,nu); B_d(5,1) = dt/tau; B_d(6,2) = dt/tau;
C_d = diag([1 1 1 1 1 1]);
D_d = zeros(ns, nu);

x_init = [0.1;0.1;0;0;0;0];

% LQ gains
gain = 0.1;
w1 = 1;
w2 = 1;
w3 = 0.01;
w4 = 0.01;
Q_n = gain*diag([w1 w2 w3 w4 0 0]);
R_n = 1e-4*eye(nu);

L_k = zeros(N,nu,ns);
S = zeros(N,ns,ns);

S(end,:,:) = Q_n;

for i = N:-1:2
    L_k(i-1,:,:) = (R_n + B_d'*squeeze(S(i,:,:))*B_d)\B_d'*squeeze(S(i,:,:))*A_d;
    S(i-1,:,:) = Q_n + A_d'*squeeze(S(i,:,:))*(A_d-B_d*squeeze(L_k(i-1,:,:)));
end

% Kalman
Sigma = zeros(N, ns, ns);
K_k = zeros(N, ns, ns);

Sigma(i,:,:) = 1e-7*eye(ns);
oXi = 0.1*(B_d*B_d');
oOmega = 0.1*max(max(oXi))*eye(ns);

for j = 1:N-1
   K_k(j,:,:) = A_d*squeeze(Sigma(j,:,:))*C_d'/(C_d*squeeze(Sigma(j,:,:))*C_d' + oOmega);
   Sigma(j+1,:,:) = oXi + (A_d-squeeze(K_k(j,:,:))*C_d)*squeeze(Sigma(j,:,:))*A_d';
end

delay = 10;