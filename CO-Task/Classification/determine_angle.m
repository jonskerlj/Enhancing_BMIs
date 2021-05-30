% Function that converts target angle into labels from 1 to 8
function class = determine_angle(td)
   switch td.target_direction
       case 0
           class = 1;
       case pi/4
           class = 2;
       case pi/2
           class = 3;
       case 3*pi/4
           class = 4;
       case pi
           class = 5;
       case -3*pi/4
           class = 6;
       case -pi/2
           class = 7;
       case -pi/4
           class = 8;     
   end  
end