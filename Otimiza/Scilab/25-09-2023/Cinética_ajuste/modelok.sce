// Modelo matemático --------------------------------------------------
function Cout = modelok(T,C0,t0,tscan)
global T
clear yy Cout
    tt  = tscan
    x0 = C0
    yy = ode(x0,t0,tt,fkinetic01)
    Cout = yy'
endfunction


// Modelo matemático --------------------------------------------------
// reação 2A -> B -> C
// reator batelada.
function ydot = fkinetic01(t,y)
global T kc
    CA = y(1)
    CB = y(2)
    CC = y(3)
    k1 = kc(1) //kc0(1)*exp(-Ea(1)/T)
    k2 = kc(2) //kc0(2)*exp(-Ea(2)/T)
    dCAdt = -k1*CA^2
    dCBdt =  k1*CA^2 - k2*CB
    dCCdt =  k2*CB
    ydot = [dCAdt dCBdt dCCdt]'
endfunction
//---------------------------------------------------------------------
