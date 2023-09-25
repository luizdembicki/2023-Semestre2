//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// UFPR
// Engenharia Química
// Disciplina de Otimização de Processos (TQ094) 
// Prof. Marcos L. Corazza

// ---------------------------------------------------------------------
// Método de BFGS: multidimensional sem restrições
// ---------------------------------------------------------------------


function [xopt,fopt,i,gx] = optBFGS(func,x0,niter,tol)

disp('---------------------------------------------------- ')
disp(' Metodo de BFGS multidimensional sem restrições ')
disp('---------------------------------------------------- ')
// Arquivo para gravação dos resultados parciais.
fd1 = mopen('saidasBFGS.txt','wt')

n = length(x0)

//Cálculo das derivadas numericas de 1a e 2a ordem da função em x.
x   = x0
fx  = func(x)
dfx = gradf(n,func,x)
gx  = dfx
Hk = eye(n,n) //inicia com a Matriz identidade.

//Imprime os resultados de cada iteração.
z = fprintr(n,0,x,fx,gx)

i=0
while (norm(gx)>tol)
  
  i = i + 1
  
  x0   = x
  gx0 = gx

  // Resolve eq. (5.15): H.s = -g
  s = descent_dir(Hk,gx)
  
  // Atualiza x | f(x+dx)< f(x) (Usa condição de Armijo)
  [x,fx,itback] = backtracking(func,x,fx,gx,s) 
  dfx = gradf(n,func,x)
  gx = dfx
  
  yk = gx - gx0  //delta g(k)
  
  dk = x - x0    //delta x(k)
  
  H1 = yk*yk'/(dk'*yk)
  H2 = Hk*dk*(Hk*dk)'/(dk'*Hk*dk)
  Hk = Hk + H1 - H2 
  //aux1 = Hk*dk*dk'*Hk/dotprod(dk,Hk*dk)
  //dHk  = (yk*yk')/dotprod(dk,yk) - aux1   // delta H(k)
  //Hk = Hk + dHk   // Atualiza a Hessiana.

  //Imprime os resultados de cada iteração.
  z = fprintr(n,i,x,fx,gx)

end
xopt = x
fopt = func(x)
mclose(fd1)
endfunction
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// --------------------------------------------------------------------
// Subrotinas auxiliares.
//
// ---------------------------------------------------------------------
// Cálculo do gradiente de uma função f(x).
// Método de Diferenças Centrais.
// ---------------------------------------------------------------------
function df = gradf(n,func,x)
   eps = 1e-5    // h tamanho da perturbação 
   xx = x           // guarda x em xx
   for i=1:n
     //if x(i)>1, h = abs(x(i))*eps, end
     //if x(i)<1, h = eps, end
     h = abs(x(i))*eps
     if h<eps, h=eps,end
     x(i) = xx(i) + h      // para trás
     f1 = func(x)          // calc. da função para frente.
     x(i) = xx(i) - h
     f2 = func(x)        // calc. da função para trás
     df(i) = (f1 - f2)/(2*h)  // aplica fórmula de dif. centrais.
     x(i) = xx(i)
   end
endfunction
// ---------------------------------------------------------------------



// ---------------------------------------------------------------------
// Calculo de direção descendente.
// ---------------------------------------------------------------------
function [d,gx] = descent_dir(Bk,gx)
  d = linsolve(Bk,gx)
endfunction
// ---------------------------------------------------------------------


// ---------------------------------------------------------------------
// line search by backtracking until Armijo condition
// ---------------------------------------------------------------------
function [xnew,fnew,itback]=backtracking(f,x,fx,gx,d) 
tau=0.3;
bet=0.0001;
alphainit=1;
alpha=alphainit;xnew=x+alpha*d;
fnew=f(xnew);
itback=1;
while(fnew>fx+bet*alpha*dotprod(gx,d))
  alpha=tau*alpha;
  xnew=x+alpha*d;
  fnew=f(xnew);
  itback=itback+1;
end
endfunction
//---------------------------------------------------------------------


//---------------------------------------------------------------------
// computes the dot product of x and y
//---------------------------------------------------------------------
function z=dotprod(x,y)
z=sum(x.*y);
endfunction
//---------------------------------------------------------------------


//-------------------------------------------------------------------
// Subrotina para impressão de resultados parciais.
//-------------------------------------------------------------------
function ff = fprintr(n,iter,xx,fxx,dfxx)
  mfprintf(fd1,'%s \n',' ')
  mfprintf(fd1,'%s %i\n','iter:',iter)
  for i=1:n
    mfprintf(fd1,'%s %i %s %f\n','x(',i,')=',xx(i))
  end
  mfprintf(fd1,'%s %f\n','f(x):',fxx)
  for i=1:n
    mfprintf(fd1,'%s %i %s %f\n','gradf(',i,')=',dfxx(i))
  end
  mfprintf(fd1,'%s \n','-------------------')
  mfprintf(fd1,'%s \n',' ')
  ff = 'done'
  //mclose(fd1)
endfunction
//---------------------------------------------------------------------


