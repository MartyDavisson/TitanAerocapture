function f = orbdynTitan( t,X )
mu = 8.977972416e+12;
ae = 2574.7e3;

x = X(1);
y = X(2);
z = X(3);
r = sqrt(x^2 + y^2 + z^2);

f = [ X(4);
      X(5);
      X(6);
      mu*((-1)*x/(r^3));	  
      mu*((-1)*y/(r^3)); 
      mu*((-1)*z/(r^3))];
end