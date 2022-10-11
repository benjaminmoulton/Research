clc;

syms x y z u v w p b L c_t c_r t_t t_r s_k s_p s_L s_T;
% u v w are x y z cg locations

x_up =  (1/4 - s_L) * ( (c_t-c_r)/b * y + c_r ) - s_k - s_p;
x_lo = -(3/4 - s_T) * ( (c_t-c_r)/b * y + c_r ) + s_k + s_p;
z_up =  1/2 * ( (t_t*c_t-t_r*c_r)/b * y + t_r*c_r ) - s_k - s_p;
z_lo = -1/2 * ( (t_t*c_t-t_r*c_r)/b * y + t_r*c_r ) + s_k + s_p;
y_up = b;
y_lo = 0;

m = simplify(int(int(int(p,x,x_lo,x_up),z,z_lo,z_up),y,y_lo,y_up));
% disp(m);
m_latex = latex(m);
% disp(m_latex);