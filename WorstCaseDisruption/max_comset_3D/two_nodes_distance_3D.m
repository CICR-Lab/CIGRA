function dist=two_nodes_distance_3D(rs,nodeA,nodeB)
R=6378.137;
% x1=R*cosd(rs.node_data(nodeA,7))*sind(rs.node_data(nodeA,6));y1=R*cosd(rs.node_data(nodeA,7))*cosd(rs.node_data(nodeA,6));z1=R*sind(rs.node_data(nodeA,7));
% x2=R*cosd(rs.node_data(nodeB,7))*sind(rs.node_data(nodeB,6));y2=R*cosd(rs.node_data(nodeB,7))*cosd(rs.node_data(nodeB,6));z2=R*sind(rs.node_data(nodeB,7));
x1=rs.node(nodeA,4);y1=rs.node(nodeA,5);z1=rs.node(nodeA,6);
x2=rs.node(nodeB,4);y2=rs.node(nodeB,5);z2=rs.node(nodeB,6);
dist=sqrt((x1-x2)^2+(y1-y2)^2+(z1-z2)^2);
dist=2*R*asin(dist/(2*R));
