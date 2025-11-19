function inter_node=solve_two_curves_equations_3D(curve1,curve2)
R=6378.137;
a1=curve1(1);a2=curve1(2);a3=curve1(3);k1=curve1(4);
b1=curve2(1);b2=curve2(2);b3=curve2(3);k2=curve2(4);
c1=(a1*b3-a3*b1)/(b2*a3-a2*b3);d1=(k2*a3-k1*b3)/(b2*a3-a2*b3);
if a3==0 && b3==0
    inter_node=[0 0 R;0 0 -R];
else 
  if a3~=0
     e1=(-a1-a2*c1)/a3;f1=(-a2*d1+k1)/a3;
  else
    e1=(-b1-b2*c1)/b3;f1=(-b2*d1+k2)/b3; 
  end
  if (c1*d1+e1*f1)^2-(1+c1^2+e1^2)*(d1^2+f1^2-R^2)>=0
      inter_node=zeros(2,3);
      x=(-(c1*d1+e1*f1)+sqrt((c1*d1+e1*f1)^2-(1+c1^2+e1^2)*(d1^2+f1^2-R^2)))/(1+c1^2+e1^2);
      y=c1*x+d1;z=e1*x+f1;inter_node(1,:)=[x,y,z];
      x=(-(c1*d1+e1*f1)-sqrt((c1*d1+e1*f1)^2-(1+c1^2+e1^2)*(d1^2+f1^2-R^2)))/(1+c1^2+e1^2);
      y=c1*x+d1;z=e1*x+f1;inter_node(2,:)=[x,y,z];
  else
      inter_node=[];
  end
end
