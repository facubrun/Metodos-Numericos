function parte2
 t = load("tiempos").tiempos;
 t22 = load("tiempos2").tiempos2;
 tiempocorregido = sort([-1,t]);
 t2 = [-1,t22];
  #MEDICIONES DE y0:
  y0=[0.9142,0.8655,1.0861,1.0832];

  hold on;
  title("Lsode para la condicion inicial y0=0.9142");
  xlabel("tiempo");
  ylabel("y");
  plot(tiempocorregido,lsodeP2_1(y0,tiempocorregido),"r")
  hold off;

 #PARTE 2
  function p2_1 = lsodeP2_1(y_0,tiempoc)
    function ecuacion= g(y,tiemp2)
      ecuacion = [(y^2)*(1-3*tiemp2)-3*y];
    endfunction
    L1 = lsode("g",y0(1),tiempoc); #para cada condicion inicial y0
    p2_1 = L1;
  endfunction

  function p2_2 = lsodeP2_2(y_0,tiempoc)
    function ecuacion= g(y,tiemp2)
      ecuacion = [(y^2)*(1-3*tiemp2)-3*y];
    endfunction
    L2 = lsode("g",y0(2),tiempoc);
    p2_2 = L2;
  endfunction


  function p2_3 = lsodeP22_3(y_0,tiempoc)
     function ecuacion2 = l(u,ti) #cambio de variable u = y^-1
       ecuacion2= [3*u+3*ti-1];
     endfunction
     y_0 = 1/y0(3);
     L3 = lsode("l",y_0,tiempoc);
     for i=1:1:152
       L3(i)=(1/L3(i));
     endfor
     p2_3 = L3;
  endfunction

  function p2_4 = lsodeP22_4(y_0,tiempoc)
     function ecuacion2 = l(u,ti) #cambio de variable u = y^-1
       ecuacion2= [3*u+3*ti-1];
     endfunction
     y_0 = 1/y0(4);
     L4 = lsode("l",y_0,tiempoc);
     for i=1:1:152
       L4(i)=(1/L4(i));
     endfor
     p2_4 = L4;
  endfunction

endfunction
