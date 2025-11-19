function dist=edge_center_distance_3D(rs,edgeA,center)
dist=Inf;R=6378.137;
n1=rs.edge_data(edgeA,2);x1=rs.node(n1,4);y1=rs.node(n1,5);z1=rs.node(n1,6);
n2=rs.edge_data(edgeA,3);x2=rs.node(n2,4);y2=rs.node(n2,5);z2=rs.node(n2,6);
% x1=R*cos(rs.node_data(n1,7))*sin(rs.node_data(n1,6));y1=R*cos(rs.node_data(n1,7))*cos(rs.node_data(n1,6));z1=R*sin(rs.node_data(n1,7));
% x2=R*cos(rs.node_data(n2,7))*sin(rs.node_data(n2,6));y2=R*cos(rs.node_data(n2,7))*cos(rs.node_data(n2,6));z2=R*sin(rs.node_data(n2,7));

dist12=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);dist12=2*R*asin(dist12/(2*R));

if dist12<=10^-4
   dist=sqrt((x1-center(1))^2+(y1-center(2))^2+(z1-center(3))^2);dist=2*R*asin(dist/(2*R));
else
   dist1=sqrt((x1-center(1))^2+(y1-center(2))^2+(z1-center(3))^2);dist1=2*R*asin(dist1/(2*R));
   dist2=sqrt((x2-center(1))^2+(y2-center(2))^2+(z2-center(3))^2);dist2=2*R*asin(dist2/(2*R)); 
   if abs(dist1)<10^-4 || abs(dist2)<=10^-4
      dist=0;
   else
      g1=y1*z2-y2*z1;g2=x2*z1-x1*z2;g3=x1*y2-x2*y1;h1=g2*center(3)-g3*center(2);h2=g3*center(1)-g1*center(3);h3=g1*center(2)-g2*center(1);
      inter_node=solve_two_curves_equations_3D([g1 g2 g3 0],[h1 h2 h3 0]);
      for vt=1:2
          inx=inter_node(vt,1);iny=inter_node(vt,2);inz=inter_node(vt,3);
          dist1in=sqrt((x1-inx)^2+(y1-iny)^2+(z1-inz)^2);dist1in=2*R*asin(dist1in/(2*R));
          dist2in=sqrt((x2-inx)^2+(y2-iny)^2+(z2-inz)^2);dist2in=2*R*asin(dist2in/(2*R));
          if abs(dist12-(dist1in+dist2in))<=0.001
             sld=sqrt((inx-center(1))^2+(iny-center(2))^2+(inz-center(3))^2);dist=2*R*asin(sld/(2*R));   
          else    
             dist=min([dist1 dist2 dist]); 
          end
      end
   end   
end


