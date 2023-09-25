
// Inicialização do arquivo.
clear
clc


// ---------------------------------------------------------------------
// Exemplo 2
function f = funcao2(x)
  f = (x(1)-2)^4+(x(1)-2)^2*x(2)^2+(x(2)+1)^2
endfunction
// ---------------------------------------------------------------------



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Exemplo de Chamada: Métodos de Levenberg-Marquard.
  exec('OptLM.sce');
  
  tol = 5.e-6;
  maxiter = 1000;
  
  x0 = [-10; 10];
  [xopt,fopt,iter] = OptLM(funcao2,x0,maxiter,tol);

  disp(xopt)
  disp(fopt)
  disp(iter)

 
 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
