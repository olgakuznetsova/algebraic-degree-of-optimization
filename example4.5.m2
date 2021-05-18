-- we compute the closest point x on a fixed ellipse to a data point u wrt 3-norm 
-- since we are taking 3rd power of 3-norm, we need to solve 2 different optimization problems
restart
KK = QQ;
R = KK[x,y];
I = ideal(x^2+4*y^2-1); -- ellipse
jacI = transpose jacobian I;
for i in 1..2 do u_i = random(KK) -- this is the random data point
u_1 = -6/10 -- special data point
u_2 = 6/10
-- when we optimize the 3-norm sqrt_3(|x-u1|^3+|y-u2|^3) with (x,y) on the ellipse and (u1,u2) given,
-- and we want to do it algebraically, then we need to split our problem into two different subproblems,
-- with the two objective functions
-- (x-u1)^3+(y-u2)^3 and (x-u1)^3-(y-u2)^3
J = ideal((x-u_1)^3+(y-u_2)^3); -- case 1
J = ideal((x-u_1)^3-(y-u_2)^3); -- case 2
jacJ = transpose jacobian J;
Icrit = I+minors(codim(I)+1,jacI||jacJ);
degree Icrit, codim Icrit
-- this is a zero-dimensional variety, and the cardinality is the 3-degree of the variety.
S = CC[x,y];
gg = first entries gens Icrit;
F = apply(gg, i-> sub(i,S));
needsPackage "Bertini"
sol = bertiniZeroDimSolve F;
coordsol = apply(sol, i-> coordinates i);
c1 = coordsol#1 -- adjust #i to correct i on each case
c2 = coordsol#3

-- the varieties with the following equations are tangent to the ellipse at (c1,c2) 
toString((x-u_1)^3+(y-u_2)^3-(c1#0-u_1)^3-(c1#1-u_2)^3) -- case 1
toString((x-u_1)^3+(y-u_2)^3-(c2#0-u_1)^3-(c2#1-u_2)^3)

toString(-(x-u_1)^3+(y-u_2)^3+(c1#0-u_1)^3-(c1#1-u_2)^3) -- case 2
toString(-(x-u_1)^3+(y-u_2)^3+(c2#0-u_1)^3-(c2#1-u_2)^3)
