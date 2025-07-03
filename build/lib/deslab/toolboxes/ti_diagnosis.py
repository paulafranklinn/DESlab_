# Biblioteca de funções de Diagnosticabilidade em TIA
# Autor: Christiano Henrique Rezende
# Orientador: João Carlos Basilio
# Coorientador: Gustavo Sousa Viana

from deslab import *
import portion as P
from deslab.toolboxes.ti_functions import *

syms('N Y')

## Renomeia os eventos observáveis no fim de um Caminho Detectável
## e os estados de Glt com os eventos que levam a ele.
## ret = 'GRF' : retorna Glt com os estados renomeados (Renomeado Final)
## ret = 'GRI' : retorna Glt com as transições renomeadas (Renomeado Intermediario)
def rename_glt(Glt,ret = 'GRF'):
    
    G,Mu = Glt
    Xr = list.copy(list(G.X))
    X0r = list.copy(list(G.X0))
    Er = list.copy(list(G.Sigma))
    Eor = list.copy(list(G.Sigobs))
    Tr = list.copy(list(transitions(G)))
    Muri = dict.copy(Mu)
    Murf = dict.copy(Mu)

    ## Calcula os DP de cada estado
    for x in G.X:
        for p in DP(Glt,x):
            ## Se o tamanho do DP for maior que 1, é certeza que 
            ## existe evento não observável
            if len(p) > 1:
                # Primeiro cria o evento novo com a concatenação
                # de todos os eventos não observáveis
                evr = str()
                for t in p:
                    evr += t[1]

                if ret == 'GRI':
                    Eor.append(evr)
                    Er.append(evr)

                ## Caso a transição observável seja um self loop, precisa criar
                ## um estado duplicado e fazer a transição renomeada para ele
                if(p[-1][0] == p[-1][2]):
                    x_ext = str(p[-1][2])+"e"
##                    print("Duplicando estado: " + str(x_ext))
                    Xr.append(x_ext)

                    Tr.remove(p[-1])
                    
                    tr2 = (x_ext,p[-1][1],x_ext)
                    Tr.append(tr2)
##                    print("Adicionando transição " + str(tr1) + " com tempo " + P.to_string(Mu[p[-1]]))
##                    print("Adicionando transição " + str(tr2) + " com tempo " + P.to_string(Mu[p[-1]]))

                    if ret == 'GRI':
                        tr1 = (p[-1][0],evr,x_ext)
                        Tr.append(tr1)
                        
                        Muri[tr1] = Mu[p[-1]]
                        Muri[tr2] = Mu[p[-1]]
                        Muri.pop(p[-1])
                    elif ret == 'GRF':
                        xr = (x_ext,evr,Mu[p[-1]])
                        Xr.append(xr)
                        tr1 = (p[-1][0],p[-1][1],xr)
                        Tr.append(tr1)
                        
                        Murf[tr1] = Mu[p[-1]]
                        Murf[tr2] = Mu[p[-1]]
                        Murf.pop(p[-1])

                ## Se não for um self loop
                else:
                    ## Renomeando transições que levam ao último estado
                    ## e renomeia o estado caso seja GRF
                    if ret == 'GRF':
                        xr = (p[-1][2],evr,Mu[p[-1]])
##                        print("Renamed State: " + str(xr))
                        tr_to = (p[-1][0],p[-1][1],xr)
##                        print("Renamed Transition To: " + str(tr_to))

                        ## Renomeando transições que partem do último estado renomeado
                        for t in G.transitions():
                            if t[0]==p[-1][2]:
##                                print("Comparando estado origem de " + str(t) + " com estado destino da transição " + str(p[-1]))
                                ## Se for uma transição de self loop
                                if t[0] == t[2]:
    ##                                x_ext = t[2]+"'"
                                    tr_from = (xr,t[1],t[2])
    ##                                print("Criando estado: " + str(x_ext))
    ##                                Xr.append(x_ext)

                                    Tr.append(tr_from)
##                                    Tr.remove(t)
    ##                                Tr.append((x_ext,t[1],t[2]))

                                    Murf[tr_from] = Mu[t]
##                                    Murf.pop(t)
    ##                                Mur[(x_ext,t[1],t[2])] = Mu[t]
                                else:
##                                    print("Que transição é essa? " + str(t))
                                    tr_from = (xr,t[1],t[2])
                                    if t in Tr:
                                        Tr.remove(t)

                                    if p[-1][2] in Xr:
##                                        print("Removendo estado " + str(p[-1][2]))
                                        Xr.remove(p[-1][2])
                                    Tr.append(tr_from)

                                    Murf[tr_from] = Mu[t]
                                    if t in Murf:
                                        Murf.pop(t)
##                                    print("Renamed Transition From: " + str(tr_from))
                                

    ##                    tr_from = [((xo,evr,Mu[(xo,e,xd)]),e,xd) for (xo,e,xd) in  if xo==p[-1][2]]

                        Xr.append(xr)
                    ## Apenas renomeia as transições que levam ao último estado
                    elif ret == 'GRI':
                        tr_to = (p[-1][0],evr,p[-1][2])

##                    print("Removendo transição " + str(p[-1]))
                    if p[-1] in Tr:
                        Tr.remove(p[-1])
                    Tr.append(tr_to)

##                    for t in tr_from:
##                        Tr.remove
                    if ret == 'GRI':
                        Muri[tr_to] = Mu[p[-1]]
                        if p[-1] in Muri:
                            Muri.pop(p[-1])
                    elif ret == 'GRF':
                        Murf[tr_to] = Mu[p[-1]]
                        if p[-1] in Murf:
                            Murf.pop(p[-1])

##    print("Renamed States: " + str(Xr))
##    print("")
##    print("Renamed Transitions: " + str(Mur))
##    print("")

    if ret == 'GRI':
        glri = ac(fsa(Xr,Er,Tr,X0r,Sigobs = Eor,name="$G_{lT}^{Ri}$"))
        glri.setgraphic(style='observer')
        glt_ri = tia(glri,Muri)
        return glt_ri
    
    elif ret == 'GRF':
##        print("Estados de GRF: " + str(Xr))
##        print("Estados iniciais de GRF: " + str(X0r))
        glrf = ac(fsa(Xr,Er,Tr,X0r,Sigobs = Eor,name="$G_{lT}^R$"))
        glrf.setgraphic(style='observer')
        glt_rf = tia(glrf,Murf)
        return glt_rf

## Descompacta os estados renomeados de Gdtr,
## renomeando os estados e criando ou substituindo as transições.
def unpack_gdtr(Gdt_renamed):

    G,Mu = Gdt_renamed
    Xr = list.copy(list(G.X))
    X0r = list.copy(list(G.X0))
    Er = list.copy(list(G.Sigma))
    Eor = list.copy(list(G.Sigobs))
    Tr = list.copy(list(transitions(G)))
##    Mur = dict.copy(Mu)
    
    xmapping = []
    Mur = dict()
    taux = list.copy(Tr)
    for xold in G.X:
        xunp = []

         

        ## Caso de estados com pelo menos um estado renomeado (xl,evr,mu), ex:
        ## (xl1, (xl2,evr2,mu2)) ou ((xl1,evr1,mu1), (xl2,evr2,mu2))
        if any(type(ele) == tuple for ele in xold):
            eunp = []
            retain_told = False
##            print("Taux copied to unpack: " + str(xold) + " :")
##            print("")
##            print(taux)
##            print("")

            for xl in xold:
                ## Tratamento caso o estado interno seja renomeado
                if type(xl) == tuple:
                    xunp.append(xl[0])
                    eunp.append((xl[1],xl[2]))
                ## Tratamento caso o estado interno NÃO seja renomeado
                else:
                    xunp.append(xl)
                    retain_told = True
            
##            xunp.sort()
            xunp = sorted(list(dict.fromkeys(xunp)), key=Py2Key)
            xunp = tuple(xunp)

##            print("X unp: " + str(xunp))
##            print("")
##            print("Transitions with " + str(xunp) + ":")
##            print("")
            x_to_xunp = [xo for (xo,sig,xd) in taux if xd==xold]
            for x in x_to_xunp:
                for sigt in eunp:
                    Tr.append((x,sigt[0],xold))
##                    print("Adicionei: " + str((x,sigt[0],xold)) + " em Tr.")
                    Er.append(sigt[0])
                    Eor.append(sigt[0])
##                    Mur[(x,sigt[0],xunp)] = sigt[1]
                    if (x,sigt[0],xold) in Mur:
                        Mur[(x,sigt[0],xold)] = Mur[(x,sigt[0],xold)] | sigt[1]
                    else:
                        Mur[(x,sigt[0],xold)] = sigt[1]
##                    print("Adicionei: " + str((x,sigt[0],xold)) + " : " + P.to_string(sigt[1])+ " em Mur.")
##                    print("")

            
            
            for t in taux:
                if t[0]==xold:
##                    Mur[(xunp,t[1],t[2])] = Mur[t]
                    if t in Mur:
                        Mur[t] = Mur[t] | Mu[t]
                    else:   
                        Mur[t] = Mu[t]
##                    print("Adicionei: " + str(t) + " : " + P.to_string(Mu[t])+ " em Mur.")
##                    Mur.pop(t)
##                    print("Removi: " + str(t) + " de Mur.")
##                    print("")
                    
                if (t[2]==xold) and (t[1] not in [sig for (sig,mu) in eunp]):
                    if retain_told:
##                        Mur[(x,t[1],xunp)] = Mur[t]
                        if t in Mur:
                            Mur[t] = Mur[t] | Mu[t]
                        else:   
                            Mur[t] = Mu[t]

##                        print("Adicionei: " + str(t) + " : " + P.to_string(Mu[t])+ " em Mur.")

##                        Mur.pop(t)
##                        print("Removi: " + str(t) + " de Mur.")
##                        print("")
                    else:
                        Tr.remove(t)
##                        print("Removi: " + str(t) + " de Tr.")
                        if t in Mur:
                            Mur.pop(t)
##                        print("Removi: " + str(t) + " de Mur.")
##                        print("")

##            print("")
##            print("--------------------------")

        ## Caso de estados renomeados únicos (xl,evr,mu)
        elif any(type(ele) == P.Interval for ele in xold):
##            t_to_xunp = [(xo,xold[1],xd) for (xo,e,xd) in G.transitions() if xd==xold]
            xunp = xold[0]
##            print("X unp: " + str(xunp))
##            print("")
##            print("Transitions with " + str(xunp) + ":")
##            print("")

            for t in taux:
                if t[2]==xold:
                    Tr.append((t[0],xold[1],t[2]))
##                    print("Adicionei: " + str((t[0],xold[1],t[2])) + " em Tr.")
                    Tr.remove(t)
##                    print("Removi: " + str(t) + " de Tr.")
                    Er.append(xold[1])
                    Eor.append(xold[1])
                    Mur[(t[0],xold[1],xunp)] = xold[2]
##                    print("Adicionei: " + str((t[0],xold[1],xunp)) + " : " + P.to_string(xold[2])+ " em Mur.")
##                    Mur.pop(t)
##                    print("Removi: " + str(t) + " de Mur.")
##                    print("")
                if t[0]==xold:
##                    Mur[(xunp,t[1],t[2])] = Mur[t]
                    if t in Mur:
                        Mur[t] = Mur[t] | Mu[t]
                    else:   
                        Mur[t] = Mu[t]
##                    print("Adicionei: " + str(t) + " : " + P.to_string(Mu[t])+ " em Mur.")
##                    Mur.pop(t)
##                    print("Removi: " + str(t) + " de Mur.")
##                    print("")
##            print("Transições que chegam em " + str(xold) + " : " + str(t_to_xunp))
            
##            print("")
##            print("--------------------------")
        ## Caso NÃO seja um estado renomeado
        else:
            xunp = xold
            for t in taux:
                if t[2]==xold:
                    if t in Mur:
                        Mur[t] = Mur[t] | Mu[t]
                    else:   
                        Mur[t] = Mu[t]
        
##        print("X_old: " + str(xold) + " || X_unpacked: " + str(xunp))
        xmapping.append((xold,xunp))

##    print(xmapping)
##    G.renamestates(xmapping)
    
    Gdun = fsa(Xr,Er,Tr,X0r,Sigobs = Eor)
    Gdun.setgraphic(style='observer')

##    print("Renomeando keys de Mur:")
##    print("")
    keyaux = list.copy(list(Mur))
    for key in keyaux:
        keyr = list.copy(list(key))
##        print("======")
##        print("Transição: " + str(keyr))
##        print("")
        renamed = False
        for xmap in xmapping:
            if xmap[0] != xmap[1]:                
                if keyr[0] == xmap[0]:
                    keyr[0] = xmap[1]
##                    print(str(xmap[0]) + " virou " + str(xmap[1]))
                    renamed = True
                if keyr[2] == xmap[0]:
                    keyr[2] = xmap[1]
##                    print(str(xmap[0]) + " virou " + str(xmap[1]))
                    renamed = True
        if renamed:
            if tuple(keyr) in Mur:
##                print("")
##                print("Jà existe a chave " + str(tuple(keyr)) + " então vamos unir os intervalos")
##                print(P.to_string(Mur[tuple(keyr)]) + " e " + P.to_string(Mur[key]) + " = " + P.to_string(Mur[tuple(keyr)] | Mur[key]))
                Mur[tuple(keyr)] = Mur[tuple(keyr)] | Mur[key]
            else:
                Mur[tuple(keyr)] = Mur[key]
##            print("Criando item: " + str(tuple(keyr)) + " : " + P.to_string(Mur[tuple(keyr)]))
##            print("Removendo: " + str(key))
##            print("")            
            Mur.pop(key)
##    print("====== \n")
##    print(Mur)
    gdt_unp = tia(Gdun.renamestates(xmapping),Mur)

    return gdt_unp


def ti_diag(Gt,failevent,ret="GD"):
    """
    Generates the automaton for fault diagnosis in TIA
    ret variable can be defined as:
            GD (Default): Returns the diagnoser GD;
            AL: Returns the label automaton with time interval AL;
            GL: Returns the product between Gt and AL.
    -------
    Example

    syms('a b c f')

    # automaton definition G_T
    Xt1 = [0, 1, 2, 3, 4, 5]
    Et1 =[a,b,c,f]
    sigobst1 = [a,c]
    X0t1 = [0]
    Xmt1 = []
    Tt1 = [(0,a,1),(1,b,2),(2,a,3),(1,f,4), (3,c,3), (4,a,5), (5,c,5)]
    
    mut1 = {(0,a,1): P.closed(1,2),
            (1,b,2): P.closed(3,4),
            (2,a,3): P.closed(2.5,4),
            (1,f,4): P.closed(0,1),
            (3,c,3): P.closed(0.5,1.5),
            (4,a,5): P.closed(2.5,4), 
            (5,c,5): P.closed(1,2)
           }
    
    G = fsa(Xt1,Et1,Tt1,X0t1,Xmt1,Sigobs=sigobst1,name="$G_{T}$")
    GT = tia(G,mut1)

    # Diagnoser
    ti_draw(GT, 'figure')
    gdt = ti_diag(GT,f)
    ti_draw(gdt,'figure')

    """
    Xl = ['N', 'Y']
    Sigl = Gt[0].Sigma
    Sigobsl = Gt[0].Sigobs
    X0l = ['N']
    Xml = []
    Tl = [('N',failevent,'Y')]
    Tl += [('N',sig,'N') for sig in Sigl if sig != failevent]
    Tl += [('Y',sig,'Y') for sig in Sigl]

    Al = fsa(Xl,Sigl,Tl,X0l,Xml,name='$A_{l_T}$',Sigobs=Sigobsl)
    mul = {t:P.closedopen(0,P.inf) for t in Tl}

    Alt = tia(Al,mul)
    if ret == "AL":
        return Alt
    
    Glt = ti_simplify(ti_product(Gt,Alt))
    Glt[0].setgraphic(style='observer')

    if ret == "GL":
        return Glt
    else:
        pi_glt = ti_proj(Glt)
        diag = ti_equi_det(pi_glt)
        return diag

## Constroi e retorna o autômato Gscc
def ti_scc(Gt,failevent):
    """
    Generates the automaton for fault diagnosis
    based in searching for Strongly Connected Components
    -------
    Example

    syms('a b u f')

    # automaton definition G_T
    Xt1 = [0, 1, 2, 3, 4, 5]
    Et1 =[a, b, u, f]
    sigobst1 = [a,b]
    X0t1 = [0]
    Xmt1 = [ ]
    Tt1 = [(0,u,1),(1,a,2),(2,b,0),(1,f,3), (3,a,4), (4,b,5), (5,f,3)]
    
    mut1 = {(0,u,1): P.closed(2,2.5),
            (1,a,2): P.closed(3,4),
            (2,b,0): P.closed(1,4),
            (1,f,3): P.closed(2,3),
            (3,a,4): P.closed(1,2),
            (4,b,5): P.closed(2,5),
            (5,f,3): P.closed(2,4)
           }
    
    G = fsa(Xt1,Et1,Tt1,X0t1,Xmt1,Sigobs=sigobst1,name="$G_{T}$")
    GT = tia(G,mut1)

    # Timed test automaton
    gscc = ti_scc(GT,f)
    ti_draw(gscc,'figure')
    """
    Glt = ti_simplify(ti_diag(Gt,failevent,'GL'))
    Gltri = rename_glt(Glt,'GRI')
    Gltrf = rename_glt(Glt)
    Gdtr = ti_equi_det(ti_proj(Gltrf))
    Gdtun = ti_simplify(unpack_gdtr(Gdtr))
    Gscct = ext_ti_product(Gdtun,Gltri)
    return Gscct


    G, Mu = Gt
    sccs = list()
    for scc in strconncomps(G):
        ln = 0
        ly = 0
        if len(scc) > 1:
            for x in scc:
                if 'N' in x:
                    ln += 1
                if 'Y' in x:
                    ly += 1
                if ln > 0 and ly > 0:
                    sccs.append(scc)
        else:
            if any(tran[0] in scc and tran[0]==tran[2] for tran in G.transitions()):
                print("Tem selfloop em " + str(scc))
                for x in scc:
                    print("x = " + str(x))
                    if 'N' in str(x):
                        ln += 1
                    if 'Y' in str(x):
                        ly += 1
                if ln > 0 and ly > 0:
                    sccs.extend(scc)
    if len(sccs) > 0:
        if ret:
            return sccs
        else:
            return True
    else:
        return False
