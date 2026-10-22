function dydt= SirGenNetworkModel(m,matrix_phi,y,vec_beta,vec_S0,vec_I0,vec_N,gamma,zeta)


%% define the present population at node i, which is also constant in time come from rest of the nodes 

Np=zeros(m,1);                                                             % initial present population at i 
Npsum=0;
for i=1:m
for k=1:m
  Npsum=Npsum+matrix_phi(k,i).*vec_N(k); % phi(k,i)*vec_N(k,1)is the total population at vertex k that travels every day to vertex i 
end
Np(i)=Npsum;
Npsum=0;
end



% initialize the set of equations for susceptible and infected and recovered ODEs at times t

dydt = zeros(3*m,1); % this is the susceptible and infected differential equations 

% define the set of ODEs (dSdt) at times t for susceptible class 
dsdtsum=0;
for i=1:m
        
for j=1:m
  for k=1:m
     dsdtsum=dsdtsum-vec_beta(j)*matrix_phi(i,j)*y(i)*matrix_phi(k,j)*y(m+k)/Np(j); 
  end
end 
     dydt(i)=dsdtsum;
 dsdtsum=0;
end 


%% define the set of ODEs (dIdt) at times t for infected class 
didtsum=0;
%didt = zeros(2*m,1); % this is the susceptible 
for i=m+1:2*m      
for j=1:m
  for k=1:m
     didtsum=didtsum+vec_beta(j)*matrix_phi(i-m,j)*y(i-m)*matrix_phi(k,j)*y(m+k)/Np(j); 
  end
end 
 dydt(i)=didtsum-gamma*y(i);
 didtsum=0;
end 

%% define the set of ODEs (dIdt) at times t for recovered class
    for i = (2*m+1):3*m
        dydt(i) = gamma*y(i);
    end
end





