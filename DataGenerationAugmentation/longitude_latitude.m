function distance=longitude_latitude(longitude1,latitude1,longitude2,latitude2)
longitude1=longitude1*3.14/180;
latitude1=latitude1*3.14/180;
longitude2=longitude2*3.14/180;
latitude2=latitude2*3.14/180;
temp=sqrt((sin((latitude1-latitude2)/2))^2+cos(latitude1)*cos(latitude2)*(sin((longitude1-longitude2)/2))^2);
distance=2*asin(temp)*6378.137;

