% %%%%%%%%%%%%%Triplex Null Model %%%%%%%%%%%%%%%%%%%
% This code generates a null model of a undirected unweighted triplex network with 
% the same average multilinks as the original multiplex network
% This is generated with a canonical mutliplex ensemble, i.e. 
% and exponential multiplex network ensemble.
%
%It takes as an input a cell array A of elements 
%A{1} adjacency matrix of the first layer of dimension N
%A{2}  adjacency matrix of the second layer of dimension N
%A{3} adjacency matrix of the third layer of dimension N
% N total number of nodes in the multiplex
% P number of desired randomized null models
%
% The output is  a cell array C of dimension P
% The element C{n} is the n-th randomized duplex networks
% C{n}{1} if the adjacency matrix of the first layer 
%          of the n-th radomized duplex network
% C{n}{2} if the adjacency matrix of the second layer 
%          of the n-th radomized duplex network
%
% C{n}{3} if the adjacency matrix of the third layer 
%          of the n-th radomized duplex network
% This code can be redistributed and/or modified
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version.
%  
% This program is distributed ny the authors in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
%  
% If you use this code please cite 
%
% [1] G. Bianconi
% Statistical mechanics of multiplex networks: Entropy and overlap." 
% Physical Review E 87, no. 6 (2013): 062806.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C] = Triplex_Null_Model(A,N,P)

A{1}=A{1}>0;
A{2}=A{2}>0;
A{3}=A{3}>0;


k100=sum(A{1}.*(1-A{2}).*(1-A{3}))';
k010=sum(A{2}.*(1-A{1}).*(1-A{3}))';
k110=sum(A{1}.*A{2}.*(1-A{3}))';
k101=sum(A{1}.*(1-A{2}).*A{3})';
k011=sum(A{2}.*(1-A{1}).*A{3})';
k111=sum(A{1}.*A{2}.*A{3})';
k001=sum((1-A{1}).*(1-A{2}).*A{3})';
 


z100=k100/(sum(k100)+(sum(k100)==0));
z010=k010/(sum(k010)+(sum(k010)==0));
z110=k110/(sum(k110)+(sum(k110)==0));
z101=k101/(sum(k101)+(sum(k101)==0));
z011=k011/(sum(k011)+(sum(k011)==0));
z111=k111/(sum(k111)+(sum(k111)==0));
z001=k001/(sum(k001)+(sum(k001)==0));


Dold=100*(rand(N,N));
n=1;
while(n<200)
    n=n+1;
       U100=(ones(N,1)*z100' );
       U010=(ones(N,1)*z010' );
       U110=(ones(N,1)*z110' );
       U101=(ones(N,1)*z101' );
       U011=(ones(N,1)*z011' );
       U111=(ones(N,1)*z111' );
       U001=(ones(N,1)*z001' );
       D=ones(N,N) +  z100*z100'+z010*z010'+z110*z110'+  z101*z101'+z011*z011'+z111*z111'+z001*z001';
       
       w=(sum( ( U100./D - diag(diag(U100./D)) )' )' );
       z100=k100 ./ (w+(w==0));
       z100=z100.*(k100>0)+0*(k100==0);
       
       w=(sum( ( U010./D - diag(diag(U010./D)) )' )' );
       z010=k010 ./(w+(w==0)) ;
       z010=z010.*(k010>0)+0*(k010==0);
       
       w=(sum( ( U110./D - diag(diag(U110./D)) )' )' );
       z110=k110 ./(w+(w==0)) ;
       z110=z110.*(k110>0)+0*(k110==0);
       
       w=(sum( ( U101./D+ - diag(diag(U101./D)) )' )' );
       z101=k101 ./ (w+(w==0));
       z101=z101.*(k101>0)+0*(k101==0);
       
       w=(sum( ( U011./D - diag(diag(U011./D)) )' )' );
       z011=k011 ./ (w+(w==0));
       z011=z011.*(k011>0)+0*(k011==0);
       
       w=(sum( ( U111./D - diag(diag(U111./D)) )' )' );
       z111=k111 ./ (w+(w==0));
       z111=z111.*(k111>0)+0*(k111==0);
       
       w= (sum( ( U001./D - diag(diag(U001./D)) )' )' );
       z001=k001 ./(w+(w==0));
       z001=z001.*(k001>0)+0*(k001==0);
       
       Dold=D;
end


B=cell(3,1);
C=cell(P,1);
for n=1:P,
x=tril(rand(N,N));
y100=tril(ones(N,N).*(x<((z100*z100')./D)));
y110=tril((x>((z100*z100')./D)).*(x<((z100*z100'+z110*z110')./D)));
y101=tril((x>((z100*z100'+z110*z110')./D)).*(x<((z100*z100'+z110*z110'+z101*z101')./D)));
y111=tril((x>((z100*z100'+z110*z110'+z101*z101')./D)).*(x<((z100*z100'+z110*z110'+z101*z101'+z111*z111')./D)));
y010=tril((x>((z100*z100'+z110*z110'+z101*z101'+z111*z111')./D)).*(x<((z100*z100'+z110*z110'+z101*z101'+z111*z111'+z010*z010')./D)));
y011=tril((x>((z100*z100'+z110*z110'+z101*z101'+z111*z111'+z010*z010')./D)).*(x<((z100*z100'+z110*z110'+z101*z101'+z111*z111'+z010*z010'+z011*z011')./D)));
y001=tril((x>((z100*z100'+z110*z110'+z101*z101'+z111*z111'+z010*z010'+z011*z011')./D)).*(x<((z100*z100'+z110*z110'+z101*z101'+z111*z111'+z010*z010'+z011*z011'+z001*z001')./D)));


B{1}=y100+y110+y101+y111;
B{2}=y010+y110+y011+y111;
B{3}=y001+y101+y011+y111;

B{1}=B{1}+B{1}';
B{2}=B{2}+B{2}';
B{3}=B{3}+B{3}'
C{n}=B;
end



  

