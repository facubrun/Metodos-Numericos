function Obligatorio2
  t = load("tiempos").tiempos;
  t22 = load("tiempos2").tiempos2;
  tiempocorregido = sort([-1,t]);
  t2 = [-1,t22];
  #MEDICIONES DE y0:
  y0=[0.9142,0.8655,1.0861,1.0832];

## PARA PARTE 5
  datosy=[0.9526];
  for i=1:columns(t22)
    datosy(end+1)=-1/(-0.99927304*exp(3*t2(i))+t2(i));
  endfor

##  Mincuadnolineal(JR,funcionF,datosy,[1,1,1],[-1,tiempos2],100)


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
  for i=1:150
    h(i)=tiempos_z(i+1)-tiempos_z(i);
  endfor

  for i=1:150
    delta(i)= y(i+1)-y(i);
  endfor

  z(1)=0;
  z(151)=0;
  for i=1:149
    z(i+1) = ((6*delta(i+1))/(h(i+1))^2)+((6*delta(i))/(h(i))^2);
  endfor

  #Matriz=zeros(152,152)

  #diagonal inferior
##  for i=2:151
##    for j=1:150
##      Matriz(i,j) = 2/h(i-2);
##    endfor
##  endfor

  a(150)=0;
  a(151)=0;
  for i=1:149
      a(i)=2/h(i);
    endfor

  #diagonal
##  Matriz(1,1)=1;
##  Matriz(152,152)=1;
##  for i=2:151
##    for j=2:151
##      Matriz(i,j) = (4/h(i-2))+ 4/h(i-1);
##      b(i-2)= Matriz(i,j);
##    endfor
##  endfor

  b(1)=1;
  b(151)=1;
  for i=2:150
    b(i)= (4/h(i-1))+ 4/h(i);
  endfor


  #diagonal superior
##  for i=2:151
##    for j=3:152
##      Matriz(i,j) = 2/h[i-1];
##      c(i-2)=Matriz(i,j);
##    endfor
##  endfor

 c(1)=0;
 c(151)=0;
 for i=2:150
      c(i)=2/h(i);
    endfor

## ALGORITMO DE THOMAS

N = 151;

% Modificamos los primeros coeficientes
c(1) = c(1) / b(1);
z(1) = z(1) / b(1);

for n = 2:1:N
    temp = b(n) - a(n) * c(n - 1);
    if (n<N)
        c(n) = c(n) / temp;
    endif
    z(n) = (z(n) - a(n) * z(n - 1)) / temp;
endfor

% sustituimos
 x(N) = z(N);
 for n = (N - 1):-1:1
    x(n) = z(n) - c(n) * x(n + 1);
 endfor


 MatrizS = zeros(4,4);
 b = zeros(4,1);
 Matrizalfas = zeros(151,4);

 for i=1:150


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


 hold on;
 for i=1:150
   intervalo = linspace(tiempos_z(i),tiempos_z(i+1),15);
   for j=1:20
    x= intervalo(j);
   polinomio = Matrizalfas(i,1) + (Matrizalfas(i,2))*x + (Matrizalfas(i,3))*x^2 + (Matrizalfas(i,4))*x^3;
   plot(x,polinomio,"b");
   endfor
 endfor
hold off;
 parte4=Matrizalfas;
endfunction

  #PARTE 5
 Matrizg1 = p4(tiempocorregido,lsodeP2_1(y0,tiempocorregido));
  Matrizg2 = p4(tiempocorregido,lsodeP2_2(y0,tiempocorregido));
  for i=1:length(tiempocorregido)
    g1(i)=gcinco(Matrizg1,t2,tiempocorregido(i));
  endfor
  for i=1:length(tiempocorregido)
    g2(i)=gcinco(Matrizg2,t2,tiempocorregido(i));
  endfor


  function result=gcinco(Matriz,Tiempos,t)
    result=0;
    for i=1:(length(Tiempos)-1)
      if(t > Tiempos(i) && t < Tiempos(i+1))
        result = Matriz(i,1) + Matriz(i,2)*t + (Matriz(i,3))*(t^2) + (Matriz(i,4))*(t^3);
        break;
      endif
    endfor
  endfunction

  function resu=funcionF(x,T,T2)
    for i=1:length(T)
      aux(i)=x(1)*g1(i)^x(2)+g2(i)^x(3); #x(1)=alfa, x(2)=beta, x(3)= gamma.
    endfor
    resu=aux;
  endfunction

  function resultado=jacr(x,T,T2) #mal escrito a proposito para evitar conflictos con funciones built-in octave
    g1aux=g1(1);
    g2aux=g2(1);
    resaux=[g1aux^x(2),x(1)*log(g1aux)*g1aux^x(2),g2aux^x(3)];
    for i=2:columns(T)
      g1aux=g1(i);
      g2aux=g2(i);
      resaux=[resaux;[g1aux^x(2),x(1)*log(g1aux)*g1aux^x(2),log(g2aux)*g2aux^x(3)]];
    endfor
    resultado=resaux;
  endfunction


 function Resultado=Mincuad(A,Y)
   Yaux=transpose(A)*Y;
   Aaux=transpose(A)*A;
   matrizaux=Aaux;
   matrizaux(:,end+1)=Yaux;
   Resultado=PivoteoParcial(matrizaux);
   endfunction

##  function Resultado=Mincuad(A,Y)
##    [Q,R] = qr(A);
##    A=Q*R;
##    matriz_aux = [A,Y];
##    Resultado=PivoteoParcial(matriz_aux);
##  endfunction
4;
  #Alg Gauss-Newton
  function xfinal=Mincuadnolineal(Y,Xin,T,T2,niteraciones)
  X=Xin

  for i=1:niteraciones
    A=jacr(X,T,T2);
    size(A)
    Yk=transpose(Y-funcionF(X,T2,T));
    P=Mincuad(A,Yk)
    X=X+P;
  endfor
  xfinal=X;
  endfunction

  Mincuadnolineal(datosy,[1,1,1],tiempocorregido,t2,10)
endfunction
