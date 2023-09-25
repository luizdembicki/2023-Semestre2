// Otimização de Processos (TQ094)
// Prof. Dr. Marcos L. Corazza
// DEQ/UFPR

// Método de Newton.


// Exemplo 1 - Derivadas analíticas.


// Inicialização do arquivo.
clear
clc

// --------------------------------------------------------------------
// Exemplo 1
function f = f1(x)
 f = 10*x(1)^2 + 5*x(2)^2 - 3
endfunction

// Calculo do gradiente da Função 'funcao1'
function [df] = g1(x)
 df(1) = 20*x(1)
 df(2) = 10*x(2)
endfunction

// Calculo da Hessiana da Função 'funcao1'
function [Hs] = H1(n,func,x,h)
  Hs(1,1) = 20
  Hs(1,2) = 0
  Hs(2,1) = 0
  Hs(2,2) = 10
 endfunction
// --------------------------------------------------------------------




// --------------------------------------------------------------------
// Exemplo 2
function f = f2(x)
  f = (x(1)-2)^4+(x(1)-2)^2*x(2)^2+(x(2)+1)^2
endfunction

// Calculo do gradiente da Função 'funcao1'
function [df] = g2(x)
 df(1) = 4*(x(1)-2)^3+2*(x(1)-2)*x(2)^2
 df(2) = 2*(x(1)-2)^2*x(2)+2*(x(2)+1)
endfunction

// Calculo da Hessiana da Função 'funcao1'
function [Hs] = H2(n,func,x,h)
  Hs(1,1) = 12*(x(1)-2)^2+2*x(2)^2
  Hs(1,2) = 4*(x(1)-2)*x(2)
  Hs(2,1) = Hs(1,2)
  Hs(2,2) = 2*(x(1)-2)^2+2
 endfunction
// --------------------------------------------------------------------



// --------------------------------------------------------------------
// Exemplo 3 OBS: Possui + de um pt ótimo.
function f = f3(x)
  f = 5*x(1)^2 + 3*(x(1)*x(2)-2)^2 + 3*x(3)^2 + 5
endfunction

// Calculo do gradiente da Função 'funcao3'
function [df] = g3(x)
 df(1) = 10*x(1) + 6*x(2)*(x(1)*x(2)-2)
 df(2) = 6*x(1)*(x(1)*x(2)-2)
 df(3) = 6*x(3)
 //disp(df)
endfunction

// Calculo da Hessiana da Função 'funcao3'
function [Hs] = H3(n,func,x,h)
  Hs(1,1) = 10+6*x(2)^2
  Hs(1,2) = 12*x(1)*x(2)-12
  Hs(2,1) = Hs(1,2)
  Hs(1,3) = 0
  Hs(3,1) = Hs(1,3)
  Hs(2,2) = 6*x(1)^2
  Hs(2,3) = 0
  Hs(3,2) = Hs(2,3)
  Hs(3,3) = 6
  //disp(Hs)
 endfunction
// ---------------------------------------------------------------------




// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Exemplo de Chamada: Métodos de Newton Multidimensional.
// Carrega arquivo com a subrotina genérica do método de Newton.
  exec('/home/ldarch/Documentos/uni/2023-Semestre2/Otimiza/Scilab/25-09-2023/Metodo_Newton_Exemplos/optNewtonMulti.sce')
  
  
  tol = 5.e-6    //define critério de convergência.
  maxiter = 100  // número máximo de iterações.
  
  //x0 = [-5; 2];
  //[xopt,fopt,iter] = OptNewtonMulti(f2,g2,H2,x0,maxiter,tol);

   x0 = [10; 20];
  [xopt,fopt,iter] = OptNewtonMulti(f1,g1,H1,x0,maxiter,tol);

  disp('Solucao')
  disp(xopt)
  disp(fopt)
  disp(iter)
  disp(g1(xopt))
 
 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
