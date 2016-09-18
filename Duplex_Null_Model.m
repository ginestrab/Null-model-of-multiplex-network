% %%%%%%%%%%%%%Duplex Null Model %%%%%%%%%%%%%%%%%%%
% This code generates a null model of a undirected unweighted duplex network with 
% the same average multilinks as the original multiplex network
% This is generated with a canonical mutliplex ensemble, i.e. 
% and exponential multiplex network ensemble.
%
%It takes as an input a cell array A of elements 
%A{1} adjacency matrix of the first layer of dimension N
%A{2}  adjacency matrix of the second layer of dimension N
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

function [C] = Duplex_Null_Model(A,N,P)

A{1}=A{1}>0;
A{2}=A{2}>0;

k10=sum(A{1}.*(1-A{2}))';
k01=sum(A{2}.*(1-A{1}))';
k11=sum(A{1}.*A{2})';

 



tol=10^(-3);

z10=k10;
z01=k01;
z11=k11;



Dold=100*(rand(N,N));
n=1;
while(n<200)
    n=n+1;
       U10=(ones(N,1)*z10' );
       U01=(ones(N,1)*z01' );
       U11=(ones(N,1)*z11' );
       D=ones(N,N) +  z10*z10'+z01*z01'+z11*z11';
       
       z10=k10 ./ (sum( ( U10./D - diag(diag(U10./D)) )' )' );
       z10=(z10.*(k10>0)+0*(k10==0));
       
       z01=k01 ./ (sum( ( U01./D - diag(diag(U01./D)) )' )' );
       z01=(z01.*(k01>0)+0*(k01==0));
       
       z11=k11 ./ (sum( ( U11./D - diag(diag(U11./D)) )' )' );
       z11=(z11.*(k11>0)+0*(k11==0));
       err=sum(sum((D-Dold).^2))/N
       Dold=D;
end

p10=z10*z10'./D;
p01=z01*z01'./D;
p11=z11*z11'./D;

B=cell(2,1);
C=cell(P,1);
for n=1:P,
x=rand(N,N);

B{1}=tril(ones(N,N).*(x<((z10*z10'+z11*z11')./D)));
B{2}=tril((x<((z10*z10'+z11*z11'+z01*z01')./D)).*(x>((z10*z10')./D)));
B{1}=B{1}+B{1}';
B{2}=B{2}+B{2}';
C{n}=B;
end



  

