function parte1
  t = load("tiempos").tiempos;
  t22 = load("tiempos2").tiempos2;
  tc = sort([-1,t]);
  t2 = [-1,t22];
  #MEDICIONES DE y0:
  y0=[0.9142,0.8655,1.0861,1.0832];

  #Grafica
  hold on;
  title("Heun para la condicion inicial y0=0.9142");
  xlabel("tiempo");
  ylabel("y");
  plot(tc,heun(tc,y0(1)))
  hold off;


   # PARTE 1
  function metodo_Heun = heun(tiemp,y_a) #f=ecuacion
     function ecuacion= g(tiemp2,y)
       ecuacion=(y^2)*(1-3*tiemp2)-3*y;
     endfunction
     T =tiemp;
     Y=zeros(1,152);
     Y(1)=y_a;
     for i=1:151
       h = T(i+1) - T(i);
       p = Y(i) + h*g(T(i),Y(i));
       Y(i+1) = Y(i) + (h/2)*(g( T(i),Y(i)) + g(T(i+1),p));
     endfor
     metodo_Heun = Y;
  endfunction

   function metodo_Heun2 = heun2(tiemp,y_a) #f=ecuacion
     function ecuacion2= l(ti,u)
      ecuacion2= 3*u + 3*ti - 1;
     endfunction
     T =[-1, tiemp];
     Y=zeros(1,152);
     Y(1)=1/y_a;
     for i=1:151
       h = T(i+1) - T(i);
       p = Y(i) + h*l(T(i),Y(i));
       Y(i+1) = Y(i) + (h/2)*(l( T(i),Y(i)) + l(T(i+1),p));
     endfor
     for i=1:1:152
       Y(i)=(1/Y(i));
     endfor
     metodo_Heun2 = Y;
  endfunction
endfunction
