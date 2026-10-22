function [vec_S0,vec_I0,vec_N,y0] = InitialValueMobilitySIRNetwork(m,vec_N,v,Nv)
%v is the place where the infectious started
%vec_N = vec_N                % input of real origin's population vector
vec_S0 =vec_N;               % initially all people are susceptible which will be the input of real origin's population vector
vec_S0(v,1)=vec_S0(v,1)-Nv;            % Reduce the susceptible population at node v by the amount Nv, where the infection started
vec_I0=zeros(m,1);
vec_R0=zeros(m,1);
vec_I0(v,1)= Nv;  
y0=[vec_S0; vec_I0; vec_R0]';        % Initial # susceptible, and infected individuals 



