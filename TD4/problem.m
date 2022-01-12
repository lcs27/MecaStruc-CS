Parameters
% Generate polynomial parameters
polynomials
P=alpha;
clearvars -except P 
Parameters

for m = 5:N
    for n =5:N
        k_poly=polyint(conv(polyder(polyder(polyder(polyder(P(m,:))))),P(n,:)));
        K(m-4,n-4)=E*I*diff(polyval(k_poly,[0 L]));
        g_poly=polyint(conv(polyder(conv(polyder(P(m,:)),[-1,0,L^2])),P(n,:)));
        G(m-4,n-4)=1/2*rho*S*diff(polyval(g_poly,[0 L]));
        m_poly=polyint(conv(P(m,:),P(n,:)));
        M(m-4,n-4)=rho*S*diff(polyval(m_poly,[0 L]));
    end
end

phi_dot=0;
syms omega 
fun=det(K-phi_dot^2*G-omega^2*M);
fplot(fun)
omega=vpasolve(fun,omega,1)

