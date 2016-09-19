% %%%%%%%%%%%%%Swap Duplex Null Model %%%%%%%%%%%%%%%%%%%
% This code generates a null model of a undirected unweighted duplex network with 
% the same  multidegree sequence as the original multiplex network
% This is generated with a microcanonical multiplex ensemble
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

function [C] = Swap_Duplex_Null_Model(A,N,P)
Teq=10;
A{1}=A{1}>0;
A{2}=A{2}>0;

[I10,J10,V]=find(tril(A{1}.*(1-A{2})));
[I01,J01,V]=find(tril((1-A{1}).*(A{2})));
[I11,J11,V]=find(tril(A{1}.*A{2}));

C=cell(P,1);
nm=0;
 for nrun=1:(Teq*P),
     for ix=1:N,
         n1=ceil(numel(I01)*rand(1));
         i1=I01(n1);
         j1=J01(n1);
         n2=ceil(numel(I01)*rand(1));
         i2=I01(n2);
         j2=J01(n2);
          
         while(((1-A{1}(i1,j2))*(1-A{2}(i1,j2))*(1-A{1}(i2,j1))*(1-A{2}(i2,j1)))*(i1-i2)*(j1-j2)*(j1-i2)*(j2-i1)==0),
         n2=ceil(numel(I01)*rand(1));
         i2=I01(n2);
         j2=J01(n2);
         
         end
         A{2}(i1,j2)=1;
         A{2}(j2,i1)=1;
         A{2}(i2,j1)=1;
         A{2}(j1,i2)=1;
         A{2}(i1,j1)=0;
         A{2}(j1,i1)=0;
         A{2}(i2,j2)=0;
         A{2}(j2,i2)=0;
         J01(n1)=j2;
         J01(n2)=j1;
         
         n1=ceil(numel(I10)*rand(1));
         i1=I10(n1);
         j1=J10(n1);
         n2=ceil(numel(I10)*rand(1));
         i2=I10(n2);
         j2=J10(n2);
         
         while(((1-A{1}(i1,j2))*(1-A{2}(i1,j2))*(1-A{1}(i2,j1))*(1-A{2}(i2,j1))*(i1-i2)*(j1-j2)*(j1-i2)*(j2-i1))==0),
         
         n2=ceil(numel(I10)*rand(1));
         i2=I10(n2);
         j2=J10(n2);
         
         end
         A{1}(i1,j2)=1;
         A{1}(j2,i1)=1;
         A{1}(i2,j1)=1;
         A{1}(j1,i2)=1;
         A{1}(i1,j1)=0;
         A{1}(j1,i1)=0;
         A{1}(i2,j2)=0;
         A{1}(j2,i2)=0;
         
         J10(n1)=j2;
         J10(n2)=j1;
        
         
         
         n1=ceil(numel(I11)*rand(1));
         i1=I11(n1);
         j1=J11(n1);
         n2=ceil(numel(I11)*rand(1));
         i2=I11(n2);
         j2=J11(n2);
         
         while(((1-A{1}(i1,j2))*(1-A{2}(i1,j2))*(1-A{1}(i2,j1))*(1-A{2}(i2,j1)))*(i1-i2)*(j1-j2)*(j1-i2)*(j2-i1)==0),
         
         n2=ceil(numel(I11)*rand(1));
         i2=I11(n2);
         j2=J11(n2);
       
         
         end
         A{1}(i1,j2)=1;
         A{1}(j2,i1)=1;
         A{1}(i2,j1)=1;
         A{1}(j1,i2)=1;
         A{1}(i1,j1)=0;
         A{1}(j1,i1)=0;
         A{1}(i2,j2)=0;
         A{1}(j2,i2)=0;
         A{2}(i1,j2)=1;
         A{2}(j2,i1)=1;
         A{2}(i2,j1)=1;
         A{2}(j1,i2)=1;
         A{2}(i1,j1)=0;
         A{2}(j1,i1)=0;
         A{2}(i2,j2)=0;
         A{2}(j2,i2)=0;
         J11(n1)=j2;
         J11(n2)=j1;
     end
     
     if (mod(nrun,Teq)==0),
         nm=nm+1;
         C{nm}=A;
     end
 end





