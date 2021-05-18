-- we compute the critical points of the p-th power of the p-norm
-- on an affine ellipse
restart
KK = QQ;
n = 2;
R = KK[x_1..x_n];
I = ideal(x_1^2+4*x_2^2-1); -- ellipse in RR^2
jacI = transpose jacobian I;
singI = I + minors(codim I,jacI);
for i in 1..n do u_i = random(KK); -- this is the random data point
p = 4;
J = ideal sum(n, i-> (x_(i+1)-u_(i+1))^p); -- objective function
jacJ = transpose jacobian J;
Icrit = saturate(I+minors(codim(I)+1,jacI||jacJ),singI);
S = CC[x_1..x_n];
gg = first entries gens Icrit;
F = apply(gg, i-> sub(i,S));
needsPackage "Bertini"
sol = bertiniZeroDimSolve F;
