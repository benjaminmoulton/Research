syms Ixx Iyy Izz Ixy Iyx Ixz Izx Iyz Izy;

I = [Ixx,-Ixy,-Ixz;-Iyx,Iyy,-Iyz;-Izx,-Izy,Izz];

% disp(I);

Iinv = inv(I);

disp(Iinv);

C = Iinv .* (Ixx*Iyz*Izy - Ixx*Iyy*Izz + Ixy*Iyx*Izz + Ixy*Iyz*Izx + Ixz*Iyx*Izy + Ixz*Iyy*Izx);

disp(C);