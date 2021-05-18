restart
R = QQ[d_1,d_2,s,t_1,t_2,x_1,x_2,u_1,u_2,z]

E1 = d_1^2-s*t_1
E2 = d_2^2-4*s*t_2
G1 = x_1-t_1
G2 = x_2-t_2
G3 = u_1-d_1-t_1
G4 = u_2-d_2-t_2
H3 = (u_1-t_1)^2-s*t_1
H4 = (u_2-t_2)^2-4*s*t_2
Gz = z-1
F = t_1^2+4*t_2^2-1

-- all branches at once
I = ideal(E1,E2,G1,G2,H3,H4,Gz,F)
degree I, codim I --(32,8)
I' = eliminate({s,t_1,t_2},I)
degree I', codim I' --(24,5)
J = eliminate (z,I')
degree J, codim J --(24,4)
K = eliminate ({d_1,d_2},J)
degree K, codim K -- (6,2)

-- branch given by G3 and G4
I1 = ideal(E1,E2,G1,G2,G3,G4,Gz,F)
degree I1, codim I1 --(8,8)
I1' = eliminate({s,t_1,t_2},I1)
degree I1', codim I1' --(6,5)
J1 = eliminate (z,I1')
degree J1, codim J1 --(6,4)
K1 = eliminate ({d_1,d_2},J1)
degree K1, codim K1 -- (6,2)
K==K1

restart
R = CC[d_1,d_2,s,t_1,t_2,x_1,x_2,u_1,u_2,z] -- for F
R = CC[d_1,d_2,x_1,x_2,u_1,u_2,z] -- for F'
R = CC[d_1,d_2,x_1,x_2,u_1,u_2] -- for FJ
R = CC[x_1,x_2,u_1,u_2] -- for FK

-- F corresponds to the ideal I in the previous code
F = {d_1^2-s*t_1,d_2^2-4*s*t_2,-t_1+x_1,-t_2+x_2,-s*t_1+t_1^2-2*t_1*u_1+u_1^2,-4*s*t_2+t_2^2-2*t_2*u_2+u_2^2,z-1,t_1^2+4*t_2^2-1}

-- F' corresponds to the ideal I' in the previous code
F' = {z-1,x_1^2+4*x_2^2-1,d_2^2-x_2^2+2*x_2*u_2-u_2^2,d_1^2+4*x_2^2+2*x_1*u_1-u_1^2-1,x_1*x_2^2+16*x_2^3+8*x_1*x_2*u_1-4*x_2*u_1^2-2*x_1*x_2*u_2+x_1*u_2^2-4*x_2}

-- FJ corresponds to the ideal J in the previous code
FJ = {x_1^2+4*x_2^2-1,d_2^2-x_2^2+2*x_2*u_2-u_2^2,d_1^2+4*x_2^2+2*x_1*u_1-u_1^2-1,x_1*x_2^2+16*x_2^3+8*x_1*x_2*u_1-4*x_2*u_1^2-2*x_1*x_2*u_2+x_1*u_2^2-4*x_2}

-- FK corresponds to the ideal K in the previous code
FK = {x_1^2+4*x_2^2-1,x_1*x_2^2+16*x_2^3+8*x_1*x_2*u_1-4*x_2*u_1^2-2*x_1*x_2*u_2+x_1*u_2^2-4*x_2}

needsPackage "Bertini"

sol = bertiniPosDimSolve FK;

peek sol
-- F has 4 components of dim=2 (codim = 8) and degree 8
-- F' has 4 components of dim=2 (codim = 5) and degree 6
-- FJ has 4 components of dim=2 (codim = 4) and degree 6
-- FK is irreducible of dim=2 (codim = 2) and degree 6
