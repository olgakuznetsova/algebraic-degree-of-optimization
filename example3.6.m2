R = CC[x1,x2,x3,x4,u1,u2,u3,u4]
f = x1^3+x2^3+x3^2*x4-1
M = matrix{{(u1-x1)^2,(u2-x2)^2,(u3-x3)^2,(u4-x4)^2}}
F = {f}|(first entries gens minors(2, diff(matrix{{x1,x2,x3,x4}},f)||M));
needsPackage "Bertini"
sol = bertiniPosDimSolve F;

