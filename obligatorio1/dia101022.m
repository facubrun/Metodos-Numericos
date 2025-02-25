function PrimerObligatorioMetodos
  h = load("flujosh").flujosh;
  v = load("flujosv").flujosv;
  f = load("flujos").flujos;


  A = esquinas(rows(h), columns(v));

  b=fronteras(rows(h),columns(v));


  #eliminamos las variables conocidas
  b=substraer(A,b);
  A=delcolumnas(A);
  A(:,50)=[];

  #Y una fila de A y b probamos y da lo mismo cual
  A(1,:)=[];
  b(1)=[];
  maux=remover_ceros(A,b);
  r2=radioEspectralGS(maux(:,1:end-1), maux(end-1,:)) #radio espectral gauss seidel
  r1=radioEspectralJacobi(maux(:,1:end-1), maux(end-1,:)) #radio espectral jacobi
  puente(A,b) #calculo de flujo minimo para el puente
  sol1=PivoteoParcial(maux);
  sol2=inv(A)*b;
  for i=1:55
    sol1(i)-sol2(i) # pivoteo es igual a solucion por inversa
  endfor
  eig(inv(dame_la_diagonal(maux(:,1:end-1)))*(dame_la_diagonal(maux(:,1:end-1))-maux(:,1:end-1))) # vp de -inv(D)*(U+L)
  if max(abs(eig(inv(dame_la_diagonal(maux(:,1:end-1)))*(dame_la_diagonal(maux(:,1:end-1))-maux(:,1:end-1))))) < 1
    MJ(maux(:,1:end-1),maux(:,end),0.5*ones(55,1),1000,1e-12,1)
  else
    omega=relajamientoJ(maux(:,1:end-1),100000);
    if omega~=-1
      omega
      MJ(maux(:,1:end-1),maux(:,end),0.5*ones(55,1),1000,1e-12,omega)
    else
      printf("no se pudo relajar\n")
    endif
  endif

  eig(inv(dame_la_diagonal(maux(:,1:end-1)))*(dame_la_diagonal(maux(:,1:end-1))-maux(:,1:end-1))) # Esta mal la formula, serian los vp de inv((D+L))*(-U)
  if max(abs(eig(inv(dame_la_diagonal(maux(:,1:end-1)))*(dame_la_diagonal(maux(:,1:end-1))-maux(:,1:end-1))))) < 1
    MGS(maux(:,1:end-1),maux(:,end),zeros(55,1),1000,1e-12,1)
  else
    omega=relajamientoGS(maux(:,1:end-1),100000);
    if omega~=-1
      omega
      MGS(maux(:,1:end-1),maux(:,end),zeros(55,1),1000,1e-12,omega)
    else
      printf("no se pudo relajar\n")
    endif
  endif
 MGS(A,b,0,1000,0.000001,1)
end
#NO SE USO
#function respuesta=rango(matriz)
##  respuesta=rows(matriz)
##  for i=1:rows(matriz)
##    if matriz(i,i)~=0 #veo si el elemento de la diagonal es distinto de 0
##      for j=1:rows(matriz) #####NO DEBERÍA SER
##        if j~=i #si no estoy en la diagonal
##          matriz(j,:) -= matriz(i,:)*matriz(j,i)/matriz(i,i); #escalerizo
##        endif
##      endfor
##    elseif any(matriz(:,i)) #si el elemento de la diagonal es 0 busco otro que no lo sea
##      for k=1:rows(matriz)
##        if matriz(k,i)
##          aux=matriz(k,:);
##          matriz(k,:)=matriz(i,:);
##          matriz(i,:)=aux;
##        endif
##      endfor
##    else
##      respuesta =respuesta-1;
##    endif
##  endfor
##  endfunction

function res=remover_ceros(matriz,matrizaux) # Suma una fila no nula para no tener ceros en la diagonal de otra fila determinada
  matriz(:,end+1)=matrizaux;
  for i=1:rows(matriz)
    if matriz(i,i)==0
      for iaux=1:rows(matriz)
        if matriz(iaux,i)~=0
          matriz(i,:)+=matriz(iaux,:);
          break
         endif
       endfor
     endif
    endfor
    res=matriz;
  endfunction

# NO SE USO
##function res = escalerizar(matriz)
##  for j=1:1:rows(matriz)
##    for i=j+1:1:rows(matriz)
##      matriz(i,:) -= matriz(j,:)*matriz(i,j)/matriz(j,j);
##    endfor
##  endfor
##  res=matriz;
## endfunction

function Newmatriz=substraer(matriz,res) # Compensa en la matriz b los datos eliminados al quitar columnas
  f = load("flujos").flujos;
  for i=1:rows(f);
    res = res - matriz(:,f(i,1))*f(i,2);
  endfor
  Newmatriz=res;
endfunction


function Nuevamatrix=delcolumnas(matriz) # Elimina las columnas de los flujos que nos dan
  f = load("flujos").flujos;
  matriz(:,f(:,1))=[];
  Nuevamatrix=matriz;
endfunction

#PARTE 1:

#MATRIZ A

function A=esquinas(n,m)
  h = load("flujosh").flujosh;
  v = load("flujosv").flujosv;
  f = load("flujos").flujos;
  format long;
  A = zeros(n*m,(n-1)*m+(m-1)*n);
  for(i = 1:n) #filas
    for(j = 1:m) #columnas
      if j~=1 #horizontales
        A(j+(i-1)*m,j-1+(i-1)*(m-1))=(-1)^(i+1); #izquierda
      endif
      if j~=m
        A(j+(i-1)*m,j+(i-1)*(m-1))=(-1)^(i); #derecha
      endif
      if i~=n #verticales
        A(j+(i-1)*m,(m-1)*n+(j-1)*(n-1)+i)=(-1)^(j+1); #abajo
      endif
      if i~=1
        A(j+(i-1)*m,(m-1)*n+(j-1)*(n-1)+i-1)=(-1)^(j); #arriba
      endif
    endfor
  endfor
endfunction

#VECTOR B

function b=fronteras(n, m)
  h = load("flujosh").flujosh;
  v = load("flujosv").flujosv;
  f = load("flujos").flujos;
  b = zeros(m*n,1);
  for(j =1:m)
    for(i=1:n)
      if i==1
        b(j+m*(i-1)) += v(1,j)*(-1)^(j+1);
      elseif i==n
        b(j+m*(i-1)) += v(2,j)*(-1)^(j);
      endif
      if j==1
        b(j+m*(i-1)) += h(i,1)*(-1)^(i);
      elseif j==m
        b(j+m*(i-1)) += h(i,2)*(-1)^(i+1);
      endif
    endfor
  endfor

endfunction

#PARTE 3

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

#PARTE 4 Y 5

function r_espectralJ=radioEspectralJacobi (maux,b)
  D=diag(diag(maux)); #obtencion de la matriz diagonal
  L=D-tril(maux); #obtencion de la matriz diagonal superior L
  U=D-triu(maux);#obtencion de la matriz diagonal inferior U
  Qj = inv(D)*(L+U);
  r_espectralJ =max(abs(eigs(Qj)));
endfunction

function r_espectralGS=radioEspectralGS (maux,b)
  D=diag(diag(maux)); #obtencion de la matriz diagonal
  L=D-tril(maux); #obtencion de la matriz diagonal superior L
  U=D-triu(maux);#obtencion de la matriz diagonal inferior U
  Qgs=-inv(L-D)*U;
  r_espectralGS = max(abs(eigs(Qgs)));
endfunction

#Metodo de Jacobi
function resultado=MJ(A,b,x0,maxiter,tolerancia,w) # Si w = 1 Metodo de Jacobi, si 0<w<1 Metodo de relajacion de Jacobi

  resprevio=x0-tolerancia;
  res=x0;
  if A*x0~=b
    itera=0;
    while itera<maxiter && norm(res-resprevio) >tolerancia
      resprevio=res;
      itera = itera + 1;
      for i=1:columns(A)
        sumado=0;
        for j=1:columns(A)
          if j~=i
            sumado += A(i,j)*resprevio(j);
          endif
        endfor
        res(i)=(1-w)*resprevio(i)+w*(b(i)-sumado)/A(i,i);
      endfor
    endwhile
    if itera==maxiter
      printf('max iteraciones')
    endif
    resultado=res;
  else
    resultado=x0;
  endif
endfunction

function resultado=relajamientoJ(matriz,maxiter) # Obtiene w para que converja Jacobi
  diagonal=dame_la_diagonal(matriz);
  FPlusE=-matriz+diagonal;
  for j=rows(matriz)
    FPlusE(j,:)=FPlusE(j,:)/diagonal(j,j);
  endfor
  valprop=eig(FPlusE);
  res=-1;
  for j=1:maxiter
    w=j/maxiter;
    if max(abs(valprop*w+1-w))<1
      res=w;
      break
    endif
  endfor
  resultado=res;
endfunction

#Metodo de Gauss
function resultado=MGS(A,b,x0,maxiter,tolerancia,w) # Si w = 1 Metodo de Gauss, si 0<w<1 Metodo de relajacion de Gauss
  resprevio=x0;
  res=resprevio + tolerancia;
  if A*x0~=b
    itera=0;
    while itera<maxiter && norm(res-resprevio) >tolerancia
      resprevio=res;
      itera = itera +1;
      for i=1:columns(A)
        for j=1:columns(A)
          if j~=i
            if j<i
              sumado += A(i,j)*res(j);
            else
              sumado+=A(i,j)*resprevio(j);
              endif
          endif
        endfor
        res(i)=(1-w)*resprevio(i)+w*(b(i)-sumado)/A(i,i);
      endfor
    endwhile
    resultado=res;
  else
    resultado=x0;
  endif
endfunction


function resultado=relajamientoGS(matriz,maxiter)  # Obtiene w para que converja Gauss
  diagonal=dame_la_diagonal(matriz);
  DminusE=tril(matriz);
  Qtotal=inv(DminusE)*(triu(matriz)-diagonal);
  valprop=eig(Qtotal);
  res=-1;
  for j=1:maxiter
    w=j/maxiter;
    if max(abs(valprop*w+1-w))<1
      res=w;
      break
    endif
  endfor
  resultado=res;
endfunction


#funcion auxiliar
function resultado=dame_la_diagonal(matriz)
  resmatriz=zeros(rows(matriz),columns(matriz));
  for j =1:rows(matriz)
    resmatriz(j,j)=matriz(j,j);
  endfor
  resultado=resmatriz;
endfunction

#PARTE 6

function parte6=puente(A,b)
  puente5a6=zeros(55,1);
  puente5a6(5+4*7, end)=-1;
  puente5a6(6+5*7, end)=1;
  inversa_A=inv(A);
  for i=1:1000
    baux= b - i*puente5a6;
    x=inversa_A*baux;
    if all(mod(x,1)==0) && all(x>0)
      parte6=i;
      break
    endif
  endfor
endfunction
