
// Inicialização do arquivo.
clear
clc

// ---------------------------------------------------------------------
// Exemplo de utilização do Método BGFS no Scilab.
// ---------------------------------------------------------------------

// Exemplo 2.2
function f = f22(x)
  f = 100d0*(x(2)-x(1)^2)^2 + (1d0-x(1))^2
endfunction

// Construção da função objetivo no esquema do Scilab.
function [f, g, index]=fgeneral(x, index)
    ff = f22
    // Calcula f(x)
    f = ff(x);
    // Calcula a derivada numérica (gradienete) g(x)
    g = numderivative(ff,x)
endfunction




x0 = [-1.2  1];

[f0,g0] = fgeneral(x0)
disp(x0)
disp(f0)
disp(g0)

[fopt,xopt] = optim(fgeneral,x0)
[f,g] = fgeneral(xopt)

mprintf("x otimo:[%s]\n", strcat (string(xopt), " "));
mprintf("FO:[%s]\n", strcat (string(fopt), " "));
mprintf("Derivadas:[%s]\n", strcat (string(g), " "));

