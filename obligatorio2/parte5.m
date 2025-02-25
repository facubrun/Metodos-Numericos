
  t = load("tiempos").tiempos;
  t22 = load("tiempos2").tiempos2;
  tiempocorregido = sort([-1,t]);
  t2 = [-1,t22];

  y0=[0.9142,0.8655,1.0861,1.0832];
 

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


  #PARTE 5
  
  function result=gcinco(Matriz,Tiempos,t)
      result=0;
    for i=1:(length(Tiempos)-1)
      if(t >= Tiempos(i) && t <= Tiempos(i+1))
        result = Matriz(i,1) + Matriz(i,2)*t + (Matriz(i,3))*(t^2) + (Matriz(i,4))*(t^3);
        break;
      
    elseif (t>=Tiempos(end))
      result = Matriz(end,1) + Matriz(end,2)*t + (Matriz(end,3))*(t^2) + (Matriz(end,4))*(t^3);
      break
    endif
    endfor
  endfunction
  
  
  Matrizg = load("ys1").coef1;
  Matrizg2 = load("ys2").coef2;
  global g1=[];
  global g2=[];

  for i=1:length(t2)
    g1(i)=gcinco(Matrizg,tiempocorregido,t2(i));
  endfor
  for i=1:length(t2)
    g2(i)=gcinco(Matrizg2,tiempocorregido,t2(i));
  endfor


  

  function resu=funcionF(x,T,T2,f1,f2)
    for i=1:length(T)
      aux(i)=x(1)*f1(i)^x(2)+f2(i)^x(3); #x(1)=alfa, x(2)=beta, x(3)= gamma.
    endfor
    resu=aux;
  endfunction

  function resultado=jacr(x,T,T2,f1,f2) #mal escrito a proposito para evitar conflictos con funciones built-in octave
    g1aux=f1(1);
    g2aux=f2(1);
    resaux=[g1aux^x(2),x(1)*log(g1aux)*g1aux^x(2),g2aux^x(3)];
    for i=2:columns(T)
      g1aux=f1(i);
      g2aux=f2(i);
      col1=g1aux^x(2);
      col2=x(1)*log(g1aux)*(g1aux^x(2));
      col3=log(g2aux)*(g2aux^x(3));
      resaux=[resaux;[col1,col2,col3]];
    endfor
    resultado=resaux;
  endfunction
  3;
 function Resultado=Mincuad(A,Y)
   Yaux=transpose(A)*Y;
   Aaux=transpose(A)*A;
   matrizaux=Aaux;
   matrizaux(:,end+1)=Yaux;
   Resultado=PivoteoParcial(matrizaux);
   endfunction

  #Alg Gauss-Newton
  function xfinal=Mincuadnolineal(Y,Xin,T,T2,niteraciones,f1,f2)
  X=Xin;
  
  for i=1:niteraciones
    A=jacr(X,T,T2,f1,f2);
    Yk=transpose(Y-funcionF(X,T2,T,f1,f2));
    P=Mincuad(A,Yk);
    X=X+P;
  endfor
  xfinal=X;
  endfunction
    datosy=[0.9526];
  for i=1:columns(t22)
    datosy(end+1)=-1/(-0.99927304*exp(3*t22(i))+t22(i));
  endfor
  coeficientesfinales=Mincuadnolineal(datosy,[0.5,1,0.5],tiempocorregido,t2,10000,g1,g2)
  for i=1:length(tiempocorregido)
    ajustado(i)=coeficientesfinales(1)*(g1(i)^coeficientesfinales(2))+g2(i)^coeficientesfinales(3);
    endfor
  hold on;
  plot(t2,datosy,"r")
  plot(t2,ajustado,"b")
  title("Ajuste vs valores reales")
  xlabel("Tiempo")
  ylabel("Funciones")
  legend("Solucion teorica","Aproximacion por ajuste")
  hold off;
  printf("Parte 6\n")
  printf("0.20\n")
  coeficientesfinales(1)*gcinco(Matrizg,tiempocorregido,0.20)^coeficientesfinales(2)+gcinco(Matrizg2,tiempocorregido,0.20)^coeficientesfinales(3)
  printf("0.30\n")
  coeficientesfinales(1)*gcinco(Matrizg,tiempocorregido,0.30)^coeficientesfinales(2)+gcinco(Matrizg2,tiempocorregido,0.30)^coeficientesfinales(3)
  printf("0.35\n")
  coeficientesfinales(1)*gcinco(Matrizg,tiempocorregido,0.35)^coeficientesfinales(2)+gcinco(Matrizg2,tiempocorregido,0.35)^coeficientesfinales(3)
  printf("0.75\n")
  coeficientesfinales(1)*gcinco(Matrizg,tiempocorregido,0.75)^coeficientesfinales(2)+gcinco(Matrizg2,tiempocorregido,0.75)^coeficientesfinales(3)