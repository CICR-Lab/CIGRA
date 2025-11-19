function [sp_new]=dc_pf(sp)
%to get the dc power flow
%for example
%sp=[];
%sp_new.bus=[1 2 0 0 104 104 1;2 2 0 0 81 81 1;3 1 -103 -103 0 0 1;4 1 -82 -82 0 0 1];
%sp_new.branch=[1 2 1 29 29 1;1 3 1 75 75 1;2 3 1 46 46 1;2 4 1 30 30 1;3 4 1 18 18 1];
%result---187---128.333
%
BI=1;AD=2;MD=3;AG=4;MG=5;
%line_attribute:
%from_bus
%to_bus
%line_length
%line_susceptance
%line_actual_flow
%line_maximum_flow
LF=1;LT=2;LL=3;LS=4;AF=5;MF=6;
%initializing parameters
sp_new=sp;
nb=length(sp_new.bus(:,BI));%the number of bus
nl=length(sp_new.branch(:,BI));%the number of link
gbs=find(sp_new.bus(:,MG)~=0);
ngb=length(gbs);
ebs=find(sp_new.bus(:,MD)~=0 |sp_new.bus(:,MG)~=0);
%for isolated node, so there is no line.
if isempty(sp_new.branch)
   if sp_new.bus(1,MG)>0 && abs(sp_new.bus(1,MD))>0
       if abs(sp_new.bus(1,MD))<=sp_new.bus(1,MG)
          sp_new.bus(1,AD)=sp_new.bus(1,MD);
          sp_new.bus(1,AG)=abs(sp_new.bus(1,AD));
       else
          sp_new.bus(1,AD)=-sp_new.bus(1,MG);
          sp_new.bus(1,AG)=sp_new.bus(1,MG); 
       end 
   else
      sp_new.bus(:,AD)=0;
      sp_new.bus(:,AG)=0;
   end 
end
if ~isempty(sp_new.branch)
%let the node tag keep consistent.
   nmap=[(1:length(sp_new.bus(:,BI)))' sp_new.bus(:,BI)];
   sp_new.bus(:,BI)=nmap(:,1);
   for i=1:length(sp_new.branch(:,LF))
       sp_new.branch(i,1)=find(nmap(:,2)==sp_new.branch(i,1));
       sp_new.branch(i,2)=find(nmap(:,2)==sp_new.branch(i,2));
   end
   if ngb==0
      sp_new.bus(:,AD)=0;
      sp_new.branch(:,AF)=0;
   end
   if ngb>0
      if sum(sp_new.bus(:,AG))>=sum(abs(sp_new.bus(:,AD)))
         dl=sum(sp_new.bus(:,AG))-sum(abs(sp_new.bus(:,AD)));
         sp_new.bus(:,AG)=sp_new.bus(:,AG)*(1-dl/sum(sp_new.bus(:,AG)));  
      end
      if sum(sp_new.bus(:,AG))<sum(abs(sp_new.bus(:,AD))) && sum(sp_new.bus(:,MG))>=sum(abs(sp_new.bus(:,AD)))
         dl=sum(abs(sp_new.bus(:,AD)))-sum(sp_new.bus(:,AG));
         sp_new.bus(:,AG)=sp_new.bus(:,AG)+(sp_new.bus(:,MG)-sp_new.bus(:,AG))*dl/sum(sp_new.bus(:,MG)-sp_new.bus(:,AG));  
      end
      if sum(sp_new.bus(:,MG))<sum(abs(sp_new.bus(:,AD)))
         dl=sum(abs(sp_new.bus(:,AD)))-sum(sp_new.bus(:,MG));
         sp_new.bus(:,AG)=sp_new.bus(:,MG);
         sp_new.bus(:,AD)=sp_new.bus(:,AD)-sp_new.bus(:,AD)*dl/sum(abs(sp_new.bus(:,AD)));  
      end   
      ref=gbs(1);
      %computing the susceptance matrix
      B_matrix=sparse([sp_new.branch(:,LF);sp_new.branch(:,LT)],[sp_new.branch(:,LT);sp_new.branch(:,LF)],[-sp_new.branch(:,LS);-sp_new.branch(:,LS)]);
      for i=1:nb
          B_matrix(i,i)=-sum(B_matrix(i,:));
      end
      %F=AP,computing matrx A.
      B_matrix(ref,:)=[];
      B_matrix(:,ref)=[];
      N_matrix=sparse([(1:nl)';(1:nl)'],[sp_new.branch(:,LF);sp_new.branch(:,LT)],[sp_new.branch(:,LS);-sp_new.branch(:,LS)]);
      N_matrix(:,ref)=[];
      A_matrix=N_matrix/B_matrix;
      if ref==1
         T_matrix(:,1)=zeros(nl,1);
         T_matrix(:,2:nb)=A_matrix;
      end
      if ref==nb
         T_matrix=[A_matrix zeros(nl,1)];
      end
      if ref>1 && ref<nb
         T_matrix(:,1:ref-1)=A_matrix(:,1:ref-1);
         T_matrix(:,ref)=0;
         T_matrix(:,ref+1:nb)=A_matrix(:,ref:nb-1);
      end
      A_matrix=T_matrix(:,ebs);
      %computing the real power injection, generation-load
      rbp=sp_new.bus(:,AD)+sp_new.bus(:,AG);
      sp_new.branch(:,AF)=A_matrix*rbp(ebs);
   end
   sp_new.bus(:,BI)=nmap(:,2);
   for i=1:length(sp_new.branch(:,LF))
       sp_new.branch(i,1)=nmap(sp_new.branch(i,1),2);
       sp_new.branch(i,2)=nmap(sp_new.branch(i,2),2);
   end
end