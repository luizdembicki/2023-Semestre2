clc
clear


// Exemplo--------------------------------------------------------------
// Ajuste de parâmetros de modelos cinético
//
// Exemplo de Reação:
// 2A -->  B --> C
//   kc1    kc2

//variáveis globais
global dados npexp nc kc


//----------------------------------------------------------------------
// Sistema:  A(1) -> B(2) -> C(3)
nc = 3
comp(1) = 'A'
comp(2) = 'B'
comp(3) = 'C'

//----------------------------------------------------------------------
// Parametros do Modelo Cinético
// kc1, kc2


// Subrotinas extrnas
exec('optLM.sce',-1)    //carrega a subrotina de otimização LM.
exec('modelok.sce',-1)  //carrega subrotina com o modelo cinético.

//----------------------------------------------------------------------
//  Leitura de dados experimentais no arquivo 'dados02.txt'
// Dados do tipo Pxy
npexp = 6
dados = read('dados03.txt',npexp,3)
disp(dados)

plot(dados(:,1),dados(:,2),'ro',   dados(:,1),dados(:,3),'bo')
xlabel('t (min)')
ylabel('C [mol/L] ')
legend('CA','CB',1)

    // Resolve o modelo cinético------------------------------------
    // Dados {t,CI(0)}-> {C(t)}
       kc(1) = 2.0
       kc(2) = 0.5
       t0 = dados(1,1)
       tf = dados(npexp,1)
       C0(1) = dados(1,2)
       C0(2) = dados(1,3)
       C0(3) = 0  // Como sabemos disso? CC(0) = 0
       T = 383 // K
       t = [t0:0.1:tf]
       Csai = modelok(T,C0,t0,t)
    // -------------------------------------------------------------
plot(t,Csai(:,1),'r--',   t,Csai(:,2),'b--')

return


//----------------------------------------------------------------------
// Função de minimização.
function [ee2] = FO(p)
    global nc dados kc
   
    kc(1) = p(1)
    kc(2) = p(2)
    t     = dados(:,1)
    CAexp = dados(:,2)
    CBexp = dados(:,3)

    // Resolve o modelo cinético------------------------------------
    // Dados {t,CI(0)}-> {C(t)}
       t0 = t(1)
       C0(1) = CAexp(1)
       C0(2) = CBexp(1)
       C0(3) = 0  // Como sabemos disso? CC(0) = 0
       T = 293.15 // K
       Csai = modelok(T,C0,t0,t)
    // -------------------------------------------------------------
    
    CA = Csai(:,1)
    CB = Csai(:,2)
    CC = Csai(:,3)
    // Calcula o erro
    er1 = (CA-CAexp)^2
    er2 = (CB-CBexp)^2
    ee2 =  sum(er1) + sum(er2)
    disp(ee2)
    
endfunction



//----------------------------------------------------------------------
// Resolve a estimação dos parametros.
//----------------------------------------------------------------------


// Chute inicial nos parametros.
//kc(1) = 0.5
//kc(2) = 0.002 

disp('Ajustando ...')
tol = 1.d-5
niter = 100
p0 = kc
[par,fopt,kiter,dg] = OptLM(FO,p0,niter,tol)

disp(par)
disp(fopt)
disp(kiter)
kc(1) = par(1)
kc(2) = par(2)

    // Resolve o modelo cinético------------------------------------
    // Dados {t,CI(0)}-> {C(t)}
       t0 = dados(1,1)
       tf = dados(npexp,1)
       C0(1) = dados(1,2)
       C0(2) = dados(1,3)
       C0(3) = 0  // Como sabemos disso? CC(0) = 0
       T = 383 // K
       t = [t0:0.1:tf]
       Csai = modelok(T,C0,t0,t)
    // -------------------------------------------------------------
plot(t,Csai(:,1),'r-',   t,Csai(:,2),'b-')



// Calcula o erro médio quadrático para as variáveis usados no ajuste.
[ee2] = FO(par)
rmsd = sqrt(ee2/npexp)



// Arquivo para gravação dos resultados parciais.
fd1 = mopen('saidasAjuste.txt','wt')
mfprintf(fd1,'%s \n',' ------------------------------------------------')
mfprintf(fd1,'%s \n',' Metodo de Levemberg-Marquadt -------------------')
mfprintf(fd1,'%s \n',' ------------------------------------------------')
mfprintf(fd1,'%s \n',' ')
mfprintf(fd1,'%s %i\n','iter:',kiter)
mfprintf(fd1,'%s %f\n','kc(1)=',kc(1))
mfprintf(fd1,'%s %f\n','kc(2)=',kc(2))
mfprintf(fd1,'%s %f\n','rmsdp:',rmsd)
mfprintf(fd1,'%s %f\n','FO(par):',fopt)
for i=1:max(size(par))
    mfprintf(fd1,'%s %i %s %f\n','gradfx(',i,')=',dg(i))
end
mfprintf(fd1,'%s \n','-------------------')
mfprintf(fd1,'%s \n',' ')
mclose(fd1)


return


// Resultados: Comparação do modelo ajustado ao dados dados experimentais.
disp('Diagrama ...')

