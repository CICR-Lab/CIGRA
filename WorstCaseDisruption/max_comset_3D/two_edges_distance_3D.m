function dist=two_edges_distance_3D(rs,edgeA,edgeB)
dist=Inf;R=6378.137;
if ~isempty(intersect(rs.edge_data(edgeA,2:3),rs.edge_data(edgeB,2:3)))
    dist=0;
else
    x1=rs.node(rs.edge_data(edgeA,2),4);y1=rs.node(rs.edge_data(edgeA,2),5);z1=rs.node(rs.edge_data(edgeA,2),6);
    x2=rs.node(rs.edge_data(edgeA,3),4);y2=rs.node(rs.edge_data(edgeA,3),5);z2=rs.node(rs.edge_data(edgeA,3),6);
    x3=rs.node(rs.edge_data(edgeB,2),4);y3=rs.node(rs.edge_data(edgeB,2),5);z3=rs.node(rs.edge_data(edgeB,2),6);
    x4=rs.node(rs.edge_data(edgeB,3),4);y4=rs.node(rs.edge_data(edgeB,3),5);z4=rs.node(rs.edge_data(edgeB,3),6);
    d1=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);d2=sqrt((x3-x4)^2+(y3-y4)^2+(z3-z4)^2);
    if d1<=0.1 && d2<=0.1
       dist=sqrt((x1-x3)^2+(y1-y3)^2+(z1-z3)^2);dist=2*R*asin(dist/(2*R));    
    end
    if d1<=0.1 && d2>0.1
       dist=edge_center_distance_3D(rs,edgeB,[x1 y1 z1]);
    end
    if d1>0.1 && d2<=0.1
       dist=edge_center_distance_3D(rs,edgeA,[x3 y3 z3]);
    end
    if d1>0.1 && d2>0.1
       g1=y1*z2-y2*z1;g2=x2*z1-x1*z2;g3=x1*y2-x2*y1;k1=0; 
       h1=y3*z4-y4*z3;h2=x4*z3-x3*z4;h3=x3*y4-x4*y3;k2=0; 
       inter_node=solve_two_curves_equations_3D([g1 g2 g3 k1],[h1 h2 h3 k2]);
       for vt=1:length(inter_node(:,1))
           inx=inter_node(vt,1);iny=inter_node(vt,2);inz=inter_node(vt,3);
           dist12=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);dist12=2*R*asin(dist12/(2*R));
           dist1in=sqrt((x1-inx)^2+(y1-iny)^2+(z1-inz)^2);dist1in=2*R*asin(dist1in/(2*R));
           dist2in=sqrt((x2-inx)^2+(y2-iny)^2+(z2-inz)^2);dist2in=2*R*asin(dist2in/(2*R));
           
           dist34=sqrt((x3-x4)^2+(y3-y4)^2+(z3-z4)^2);dist34=2*R*asin(dist34/(2*R));
           dist3in=sqrt((x3-inx)^2+(y3-iny)^2+(z3-inz)^2);dist3in=2*R*asin(dist3in/(2*R));
           dist4in=sqrt((x4-inx)^2+(y4-iny)^2+(z4-inz)^2);dist4in=2*R*asin(dist4in/(2*R));
           if abs(dist12-(dist1in+dist2in))<=0.01 && abs(dist34-(dist3in+dist4in))<=0.01
              dist=0;   
           else    
              d13=edge_center_distance_3D(rs,edgeA,[x3 y3 z3]);
              d14=edge_center_distance_3D(rs,edgeA,[x4 y4 z4]);
              d21=edge_center_distance_3D(rs,edgeB,[x1 y1 z1]);
              d22=edge_center_distance_3D(rs,edgeB,[x2 y2 z2]);
              dist=min([dist d13 d14 d21 d22]);
           end  
       end
    end
end