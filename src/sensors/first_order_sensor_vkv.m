function d=first_order_sensor_vkv(a,b,N)
deltat=1/N;
x=deltat:deltat:1;
x2=x;
x2(N)=[];
x1=[0,x2];
t1=(x'*ones(1,N)-ones(N,1)*x).*b./a;
t2=(x1'*ones(1,N)-ones(N,1)*x).*b./a;
t3=(x'*ones(1,N)-ones(N,1)*x1).*b./a;
t4=(x1'*ones(1,N)-ones(N,1)*x1).*b./a;

CONST=-a./b.^2./deltat;
d=CONST.*(exp(-t1)-exp(-t2)-exp(-t3)+exp(-t4));
di=1./b+CONST.*(1-exp(-b./a.*deltat));
d=d-triu(d)+di.*eye(size(d));
