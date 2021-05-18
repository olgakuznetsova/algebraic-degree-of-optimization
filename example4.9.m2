-- script to compute the p-evolute of a plane curve
-- for p=2, this is the classic evolute, where g defines the family of normal lines
restart
R = QQ[s,t,x,y]
f = s^2+4*t^2-1
Isingf = ideal(f)+minors(1,compress transpose jacobian ideal f);
p = 4
g = det(matrix{{(x-s)^(p-1),(y-t)^(p-1)}}||diff(matrix{{s,t}},f))
h = det(diff(matrix{{s,t}},f)||diff(matrix{{s,t}},g))
Iev = eliminate({s,t}, saturate(ideal(f,g,h),Isingf));
dec = decompose Iev; -- one of the components is the p-evolute
