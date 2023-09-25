
// Inicialização do arquivo.
clear
clc

// ------------------------------------------------------------------
// Exemplo 1
function f = f1(x)
 f = 2*x(1).^2 + x(2).^2 - 3
endfunction

// --------------------------------------------------------------------
// Exemplo 2.1
function f = f21(x)
  f = x(1)^4 - 2*x(2)*x(1)^2 + x(2)^2 + x(1)^2 - 2*x(1) + 5
endfunction

// --------------------------------------------------------------------
// Exemplo 2.2
function f = f22(x)
  f = 100d0*(x(2)-x(1)^2)^2 + (1d0-x(1))^2
endfunction

// --------------------------------------------------------------------


// Exemplo 3 OBS: Tem + de um pt ótimo.
function f = f3(x)
  f = 5*x(1)^2 + 3*(x(1)*x(2)-2)^2 + 3*x(3)^2 + 5
endfunction
//---------------------------------------------------------------------


// Exemplo de Chamadas do Metodo BGFS
  
 // Todos estes Método são para problemas n-dim com derivadas numéricas.
 exec('optBFGS.sce',-1)
  

  
  tol = 1.e-5;
  maxiter = 200;
  
 
  // Exemplo f1.
   //x0 = [10.0; 10.0];
   
  // Exemplo f2.1
   //x0 = [1.0; 2.0];
   
  // Exemplo f2.2
   x0 = [-1.2d0; 1.d0];
   
  // Exemplo f3
   //x0 = [1; -2; 3]
   //x0 = [-1; 0; 0]
   //x0 = [1; -2; 1];



  [xopt,fopt,iter,df] = optBFGS(f22,x0,maxiter,tol)



disp(xopt)

mprintf('%s \n','    ')
mprintf('%s \n','    ')
mprintf('%s \n','Solução:')
mprintf('%s \n','    ')
mprintf('%s \n','xopt:')
mprintf('%f \n', xopt)
mprintf('%s \n','    ')
mprintf('%s \n',' fopt:')
mprintf('%f \n', fopt)
mprintf('%s \n','    ')
mprintf('%s \n',' gradf:')
mprintf('%f \n', df)
mprintf('%s \n','    ')
mprintf('%s \n',' iter:')
mprintf('%i \n', iter)
 
  

 
