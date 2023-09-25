//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// UFPR
// Engenharia Química
// Disciplina de Otimização de Processos (TQ094) 
// Prof. Marcos L. Corazza

// ---------------------------------------------------------------------
// Método de Newton "Puro": multidimensional sem restrições, com 
// derivadas de 1a e 2a ordem analíticas.
// ---------------------------------------------------------------------


function [xopt,fopt,k] = OptNewtonMulti(func,grad,Hessiana,x0,niter,tol)

disp('---------------------------------------------------- ')
disp(' Metodo de Newton multidimensional sem restrições '   )
disp('---------------------------------------------------- ')

x = x0
n = max(size(x))
   
clear fpv
clear xpv
k = 0
f = func(x)
g = grad(x)
disp(g)
B = Hessiana(x)
disp(B)

tt = norm(g)

while (norm(g)>tol)
 k = k + 1
 disp(k)
 dx = linsolve(B,g)
 disp(dx)
 //dx = -inv(B)*g
 //disp(dx)
 x = x + dx
 disp(x)
 f = func(x)
 g = grad(x)
 disp(g)
 B = Hessiana(x)
 if (k>niter) then
   disp('Número máximo de iterações')
   xopt = x
   fopt = f
   return
 end
end
xopt = x
fopt = f
endfunction
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++









