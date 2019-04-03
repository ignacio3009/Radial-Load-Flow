clear;clc;
format long

%Data
DATA = csvread('datos.csv');
Vb = 12.66; %kV
Sb = 10; %MVA
Zb = Vb^2/Sb; 
Nb = length(DATA); %number of branches
N = Nb + 1; %number of nodes
st=zeros(1,1000); %Stack
V = ones(1,N); %Abs(voltages) in nodes


%R and X to pu
for i=1:length(DATA)
    for j=4:5
        DATA(i,j) = DATA(i,j)/Zb;
    end
end

%P and Q to pu
for i=1:length(DATA)
    for j=6:7
        DATA(i,j) = DATA(i,j)*0.001/Sb; %Data in kW and kVAr !!!
    end
end


DVMAX=100;
tolerance = 0.0001;
count=0;
limit_count = 1e5;

% Step 1: Initializate variables
Vold = ones(1,N);
Vnew = ones(1,N);
dif = zeros(1,N);
I = zeros(1,Nb);



% Computation of branches current
while((DVMAX>tolerance && count<=limit_count))
I = zeros(1,Nb);
for j=Nb:-1:1
    flag = 0;
    for i=1:Nb
        if(DATA(j,3) == DATA(i,2)) 
            flag=1;
            break;
        end
    end
   %If it is a branch with a leaf node
   if(flag==0)
       P = DATA(j,6);
       Q = DATA(j,7);
       I(j)=(P-1i*Q)/conj(Vnew(j+1));
   else
       top=0;
       for i=2:Nb
          if(DATA(j,3)==DATA(i,2)) %Sending of i = Receiver of j
              top = top + 1;
              st(top) = DATA(i,1);
          end
       end
       while(top>0)
          temp = st(top);
          I(j) = I(j)+I(temp);
          top = top-1;
       end
       P = DATA(j,6);
       Q = DATA(j,7);
       I(j) = I(j)+(P-1i*Q)/conj(Vold(j+1));
   end
end

% Step 3:Compute new node voltages
for i=2:N
    Vold(i) = Vnew(i); %Save previous value of Vi
    Rm = DATA(i-1,4);
    Xm = DATA(i-1,5);
    Zm = Rm + 1i*Xm;
    Im = I(i-1);
    Vj = 1;
    %Look for the sending node whose receiver is node i
    for k=1:Nb
        if(DATA(k,3)==i)
            Vj = Vold(DATA(k,2));
            break;
        end
    end
    
    Vnew(i) = Vj - Im*Zm;
    count = count + 1;
end
% Step 4: Compute diference
dif = abs(Vnew - Vold);
V = abs(Vnew);
DVMAX = max(dif);
%DVMAX = tolerance; %activate for 1 iteration
end

%Print Results



