clc
clear
// Modelo matemático --------------------------------------------------
function yy = modelo(T,C0,t0,tscan)
global T
    tt  = tscan
    x0 = C0
    yy = ode(x0,t0,tt,fkinetic01)
endfunction
// Modelo matemático --------------------------------------------------
// reação 2A -> B -> C
// reator batelada.
function xdot = fkinetic01(t,x)
global T
    CA = x(1)
    CB = x(2)
    CC = x(3)
    k1 = 4000*exp(-2500/T)
    k2 = 620000*exp(-5000/T)
    dCAdt = -k1*CA^2
    dCBdt =  k1*CA^2 - k2*CB
    dCCdt =  k2*CB
    xdot = [dCAdt dCBdt dCCdt]'
endfunction
//---------------------------------------------------------------------

// Simulação.
t0 = 0              // tempo inicial
t = [t0:2:10]'    // tempo final (min)
C0(1) = 1
C0(2) = 0
C0(3) = 0
T = 293.15          // K
[Csol] = modelo(T,C0,t0,t)
// Csol(i,j)
// i: componentes
// j: pontos calculados no tempo (t0:tf)
yy = Csol';
plot(t,yy)
disp(yy)

// Arquivo para gravação dos resultados simulados.
fd1 = mopen('Report.txt','wt')
mfprintf(fd1,'%s %s %s %s\n','t(min)','CA','CB','CC')
for i=1:max(size(yy))
    mfprintf(fd1,'%f %f %f %f\n',t(i),yy(i,1),yy(i,2),yy(i,3))
end
mclose(fd1)
