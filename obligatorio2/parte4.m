function parte4
  t = load("tiempos").tiempos;
  t22 = load("tiempos2").tiempos2;
  tiempocorregido = sort([-1,t]);
  t2 = [-1,t22];

  y0=[0.9142,0.8655,1.0861,1.0832];
  plot(tiempocorregido,lsodeP2_1(y0,tiempocorregido))
  coef1=p4(tiempocorregido,lsodeP2_1(y0,tiempocorregido))
  plot(tiempocorregido,lsodeP2_2(y0,tiempocorregido))
  coef2=p4(tiempocorregido,lsodeP2_2(y0,tiempocorregido))
  title("Splines cúbicos para la condicion inicial y0=0.9142 y 0.8655")
  xlabel("tiempo")
  ylabel("y")
  legend('Lsode 1','Splines 1',"Lsode 2","Splines 2")
  hold off;
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

  function heun_p1 = parte1(tiem)
    for i=1:4
      heun(tiem,y0(i));
    endfor
  endfunction


  #PARTE 2
  function p2_1 = lsodeP2_1(y_0,tiempoc)
    function ecuacion= g(y,tiemp2)
      ecuacion = [(y^2)*(1-3*tiemp2)-3*y];
   endfunction
        L1 = lsode("g",y0(1),tiempoc); #para cada condicion inicial y0
        p2_1 = L1;
       # plot(tiempoc,p2_1,"r")
  endfunction

  function p2_2 = lsodeP2_2(y_0,tiempoc)
    function ecuacion= g(y,tiemp2)
      ecuacion = [(y^2)*(1-3*tiemp2)-3*y];
   endfunction
        L2 = lsode("g",y0(2),tiempoc);
        p2_2 = L2;
       # plot(tiempoc,p2_2,"r")
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
      # plot(tiempoc,p2_3,"r")
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
       #plot(tiempoc,p2_4,"r")
  endfunction

  #PARTE 3

  function parte3Heun = graficaHeun(x,y)  # x=valores de tiempo, y=valores de la sucesion Y(i)
    plot(x,y)
  endfunction

  function parte3Lsode = graficaLsode(w,z) # w=valores del tiempo, z = valores de la solucion dada por la funcion lsode
    plot(w,z)
  endfunction


    #PARTE 4


 function resultado=PivoteoParcial(matriz)
  filas=rows(matriz);
  for i=1:filas
    if matriz(i,i)==0
      for j=(i+1):filas
        if matriz(j,i)~=0
          auxiliar=matriz(j,:);
          matriz(j,:)=matriz(i,:);
          matriz(i,:)=auxiliar;
          break
        endif
      endfor
    endif

   for iaux=(i+1):filas
     matriz(iaux,:) -= matriz(i,:)*matriz(iaux,i)/matriz(i,i);
   endfor
 endfor
 res(filas)= matriz(filas, filas+1)/matriz(filas,filas);
 for i=filas-1:-1:1
	acum = 0;
	for p = i+1:filas
		acum = acum + (matriz(i,p)*res(p));
	endfor
  res(i)= (matriz(i, filas + 1) - acum)/matriz(i, i);
 endfor
 resultado=res;
 endfunction

function parte4=p4(tiempos_z,y)
  for i=1:151
    h(i)=tiempos_z(i+1)-tiempos_z(i);
  endfor

  for i=1:151
    delta(i)= y(i+1)-y(i);
  endfor

  z(1)=0;
  z(152)=0;
  for i=1:150
    z(i+1) = ((6*delta(i+1))/(h(i+1))^2)+((6*delta(i))/(h(i))^2);
  endfor



  a(151)=0;
  a(152)=0;

  for i=1:150
      a(i)=2/h(i);
    endfor


  b(1)=1;
  b(152)=1;
  for i=2:151
    b(i)= (4/h(i-1))+ 4/h(i);
  endfor



 c(1)=0;
 c(152)=0;
 for i=2:151
      c(i)=2/h(i);
    endfor
#Algoritmo de thomas 
N = 152;

#Modificamos los primeros coeficientes
c(1) = c(1) / b(1); 
z(1) = z(1) / b(1);

for n = 2:1:N
    temp = b(n) - a(n) * c(n - 1);
    if (n<N)
        c(n) = c(n) / temp;
    endif
    z(n) = (z(n) - a(n) * z(n - 1)) / temp;
endfor

#Sustituimos
x(N) = z(N);
for n = (N - 1):-1:1
    x(n) = z(n) - c(n) * x(n + 1);
endfor


 MatrizS = zeros(4,4);
 b = zeros(4,1);
 Matrizalfas = zeros(151,4);

for i=1:151


 MatrizS(1,:) = [1 , tiempos_z(i) , tiempos_z(i)^2 , tiempos_z(i)^3];
 MatrizS(2,:) = [1 , tiempos_z(i+1) , (tiempos_z(i+1))^2 , (tiempos_z(i+1))^3];
 MatrizS(3,:) = [0 , 1 , 2*(tiempos_z(i)) , 3*(tiempos_z(i))^2];
 MatrizS(4,:) = [0 , 1 , 2*(tiempos_z(i+1)) , 3*(tiempos_z(i+1))^2];


 b(1)= y(i);
 b(2)= y(i+1);
 b(3)= z(i);
 b(4)= z(i+1);

 matriz_parte4 = [MatrizS,b];
 alfas = PivoteoParcial(matriz_parte4);

 Matrizalfas(i,:) = alfas;

  endfor
 #falta guardar alfas de cada intervalo y graficar. endfor
 hold on;
 for i=1:151
   intervalo = linspace(tiempos_z(i),tiempos_z(i+1),25);
   for j=1:25
    x= intervalo(j);
   polinomio = Matrizalfas(i,1) + (Matrizalfas(i,2))*x + (Matrizalfas(i,3))*x^2 + (Matrizalfas(i,4))*x^3;
   plot(x,polinomio,"b");
   endfor
 endfor
 parte4=Matrizalfas;
endfunction
  

  save ys1.mat coef1
  save ys2.mat coef2
  endfunction


