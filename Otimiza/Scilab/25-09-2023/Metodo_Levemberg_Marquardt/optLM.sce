//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// UFPR
// Engenharia Química
// Disciplina de Otimização de Processos (TQ094) 
// Prof. Marcos L. Corazza

// ---------------------------------------------------------------------
// Método de Levemberg-Marquadt: multidimensional sem restrições, com 
// derivadas de 1a e 2a ordem numéricas.
// ---------------------------------------------------------------------


function [xopt,fopt,k,g] = OptLM(func,x0,niter,tol)

disp('--------------------------------------------------------------- ')
disp(' Metodo de Leverberg-Marquardt multidimensional sem restrições '   )
disp('--------------------------------------------------------------- ')

x = x0
n = max(size(x))
   
clear fpv
clear xpv
k = 0
f1 = func(x)
H = Hf1(n,func,x) //Calcula a Hessiana da função f(x)
tt = 1d30
I = eye(n,n)
B = 1000

while (tt>tol)
 k = k + 1
 g = gradf(n,func,x)
 tt = norm(g)
 H = Hf1(n,func,x) //Calcula a Hessiana da função f(x)
 Hm = (H+B*I)
 //Equação generica Hm(k).s(k) + grad(f)=0
 sk = linsolve(Hm,g)
  //ou resolve da seguinte forma:
  // Equação generica s(k) =- inv[Hm(k)].grad(f)
 //sk = -inv(Hm)*g

 alpha = 1 // ******** Ataulizar código para busca em linhas e usar p1.
 //p1 = d'*sk
 
 x = x + alpha*sk
 f2 = func(x)
 if (f2<f1) B=B/4 //,else,B=2*B,end
     f1=f2 // atualiza para a próxima iteração.
 else
     B=2*B
     x = x - sk // retorna para o x anterior.
 end
 
 if (k>=niter) then
   disp('Número máximo de iterações')
   xopt = x
   fopt = f1
   return
 end
end
xopt = x
fopt = f1
endfunction
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// ---------------------------------------------------------------------
// Cálculo do gradiente de uma função f(x).
// Método de Diferenças Centrais.
// ---------------------------------------------------------------------
function df = gradf(n,func,x)
   eps = 5e-5    // h tamanho da perturbação 
   xx = x           // guarda x em xx
   for i=1:n
     //if x(i)>1, h = abs(x(i))*eps, end
     //if x(i)<1, h = eps, end
     h = abs(xx(i))*eps
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
// Cálculo do gradiente de uma função f(x).
// Método de Diferenças Centrais.
// Opção 1.
// ---------------------------------------------------------------------
function Hs = Hf1(n,func,x)
   eps = 5e-5    // h tamanho da perturbação 
   xx = x           // guarda x em xx
   for i=1:n
       //if x(i)>1, hx = abs(x(i))*eps, end
       //if x(i)<1, hx = eps, end
       hx = abs(xx(i))*eps
       if hx<eps, hx=eps,end
       x(i) = xx(i) + hx;
       df1 = gradf(n,func,x)
       x(i) = xx(i) - hx;
       df2 = gradf(n,func,x)
       for j=1:n
         Hs(i,j) = (df1(j)-df2(j))/(2*hx)
       end
       x(i) = xx(i)
   end
endfunction
// ---------------------------------------------------------------------









