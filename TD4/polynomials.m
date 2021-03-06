A4=[2,6*L,4*3*L^2;0,6,4*3*2*L;0,0,1];
B4=[0;0;1];
alpha=[(A4\B4).'];
for n =5:N
    ls=2:n;
    A1=ls.*(ls-1).*L.^(ls-2);
	A2=ls.*(ls-1).*(ls-2).*L.^(ls-2);
    A=[A1;A2];
    for m=4:(n-1)
        alpha_mk=alpha(m-3,1:(m-1));
        Am=[];
        for l=2:n
            Am=[Am,alpha_mk*(1./(l+(2:m)+1).*L.^(l+(2:m)+1)).'];
        end
        A=[A;Am];
    end
    A=[A;zeros(1,n-1)];
    A(n-1,n-1)=1;
    B=zeros(n-1,1);
    B(n-1,1)=1;
    alpha_n=A\B;
    alpha=[alpha,zeros(n-4,1);alpha_n.'];
end
alpha=[flip(alpha,2),zeros(N-3,2)];
alpha=[zeros(4,N+1);alpha];
%figure
%hold on
%x=0:0.01:1;
%for n=4:N
    %p=alpha(n-3,:);
    %plot(x,polyval(p,x));
%end
%legend('4','5','6','7','8','9','10')
