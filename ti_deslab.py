# Biblioteca de funções de Time-Interval Discrete Event Systems

from deslab import *
import portion as P
from sys import intern
from copy import deepcopy

## Documentos dessa biblioteca estão aqui:
## https://github.com/AlexandreDecan/portion

#CLASSES AUXILIARES:
class Py2Key:

    __slots__ = ("value", "typestr")

    def __init__(self, value):
        self.value   = value
        self.typestr = intern(type(value).__name__)

    def __lt__(self, other):
        try:
            return self.value < other.value
        except TypeError:
            return self.typestr < other.typestr
#============================================================


syms('N Y')

def tia(G,mu):
    Gt = (G,mu)
    return Gt

## Adiciona os intervalos aos labels das transições e desenha
## o novo autômato Gtia (sem modificar o original)
def ti_draw(*Gt):
    if type(Gt[-1]) == str:
        style = Gt[-1]
    else:
        style = 'beamer'
    for tia in Gt:
        if type(tia) != tuple:
            break
        Gtia = tia[0]
        new_sigobs = []
        new_sig = []
        if not isitempty(Gtia):
            for trans in transitions(tia[0]):
                ev_mu = trans[1]+P.to_string(tia[1][trans])
                if trans[1] in tia[0].Sigobs:
                    new_sigobs.append(ev_mu)
                new_sig.append(ev_mu)
                Gtia = Gtia.renametransition([trans[0],(trans[1],ev_mu),trans[2]])
            Gtia = Gtia.setpar(Sigobs = new_sigobs)
            draw(Gtia,style)
        else:
            draw(fsa())
    return

## Função para determinar o Caminho Detectável de um estado em Gt
def DP(Gt,state):
    G = Gt[0]
    DP = [] 
    Eobs = list(G.Sigobs)
    E_active = list(G.Gamma(state))

    for event in E_active:
        evact = G.delta(state,event)
        if isinstance(evact,frozenset): #caso seja nao det
            #transicao com o mesmo ev indo para estados diferentes
            conj = list(evact)
            for el in conj:
                DP += [[(state,event,el)]]
        else:
            conj = evact
            DP += [[(state,event,conj)]]
        
##    print("Detectable Path in DP: " + str(DP))
    aux = DP[:]
    while len(aux) != 0:
        for p in DP:
            last_event = p[len(p) - 1][1]
            last_state = p[len(p) - 1][2]
##            print("last event of " + str(p) + " is "  + str(last_event))
##            print("last state of " + str(p) + " is "  + str(last_state))
            if p in aux:
                aux.remove(p)
            if last_event not in Eobs:
                DP.remove(p)
##                E_actives_p = list(G.Gamma(last_state))
                E_actives_p = G.Gamma(last_state)
                for ev in E_actives_p:
##                    try:
                    if type(G.delta(last_state,ev)) == frozenset:
##                        print("Identificou MAIS DE UM estado!")
                        for dest_state in G.delta(last_state,ev):
##                            print("estado destino: " + str(dest_state))
                            p_actual = p[:]
##                            print("p_atual: " + str(p_actual))
                            p_new = (last_state, ev, dest_state)
                            p_actual.append(p_new)
                            aux.append(p_actual)
                            DP.append(p_actual)
##                    except TypeError:
                    else:
##                        print("Identificou APENAS UM estado!")
                        p_actual = p[:]
                        p_new = (last_state, ev, G.delta(last_state,ev))
                        p_actual.append(p_new)
                        aux.append(p_actual)
                        DP.append(p_actual)
##    print("")
    return DP

## Faz a projeção de um Time-Interval Automaton, mas o resultado
## ainda pode ser não-deterministico.
def ti_proj(Gt):

    ## Calcula o conjunto de estados destino de eventos
    ## observáveis ou não-observáveis
    def SubsetX(G, mode = 'Obs'):
        """
        mode : if 'Obs', return Subset of states that receive a transition
               of an observable event
               if 'UnObs', return Subset of states that receive a transition
               of an unobservable event
        """
        X_new = []
        if (mode == "Obs"):
            for t in G.transitions():
                if (t[1] in E_new) and (t[2] not in X_new):
                    X_new.append(t[2])

        if (mode == "UnObs"):
            for t in G.transitions():
                if (t[1] not in E_new) and (t[2] not in X_new):
                    X_new.append(t[2])
        
        return X_new
    


    ## Função interna para gerar a projeção observável do caminho e o Intervalo gerado
    def ObservableProj(path,mu):
        #mu_temp = dict(mu)
        # Se caminho possui apenas uma transição, então ja é com um evento observavel
        # retorna a propria transição com seu intervalo original
        if len(path) == 1:            
            proj_path = path[0]
            mu_dict = {path[0] : mu.get(path[0])}
##            mu_list = [[trans, ti, io] for trans, ti, io in mu if trans == path[0]][0]
            return proj_path, mu_dict
##            return proj_path, mu_list

        else:
        # Se caminho possuir mais de uma transição, 
        # então tem que tirar a projeção observável
            initial_state = path[0][0]
            final_state = path[len(path) - 1][2]
            obs_event = path[len(path) - 1][1]
            mu_dict = dict()
##            lower_in_out = mu[path[0]].left
##            lower_in_out = [l_io for trans, ti, [l_io,u_io] in mu if trans == path[0]][0]
            
##            lower = 0
##            upper = 0
##            upper_in_out = i
            ti = P.closed(0,0)
##            print("Actual Path = " + str(path))
            for t in path:
##                t_low = mu[t][0][0]
##                t_upp = mu[t][1][0]
##                ti = [ti for trans, ti, io in mu if trans == t][0]
##                print("t = " + str(t))
                ti = interval_offset(ti,mu[t])
##                lower += ti[0]
##                upper += ti[1]
##                upper_in_out = mu[t][1][1]
##                upper_in_out = [u_io for trans, ti, [l_io,u_io] in mu if trans == t][0]

            proj_path = (initial_state, obs_event, final_state)
##            mu_list = [proj_path, [lower, upper], [lower_in_out, upper_in_out]]
            mu_dict[proj_path] = ti
            
##          return proj_path, mu_list
            return proj_path, mu_dict

##  \/ MAIN TI PROJ ================================================ \/
    G = Gt[0]
    mu_old = Gt[1]
    X0_new = list(G.X0)
    X_new = X0_new[:]
    E_new = list(G.Sigobs)
    T_new = []
    mu_new = dict()
    
    # Adiciona a X_new os estados de destino das transições rotuladas por eventos observáveis
##    print(str(X_new) + " antes da função SubsetX")
    X_new += SubsetX(G, 'Obs')
##    print(str(X_new) + " depois da função SubsetX")
    
    # Calcula o Caminho Detectável de cada estado em X com Gamma > 0
    path = []
    for x in X_new:
        if list(G.Gamma(x)):
            dp = DP(Gt,x)
            path += dp
##    print("All Paths = " + str(path))
            
    # Projeta os Caminhos Detectáveis, criando as novas transições e seus I associados
    proj = []
    for p in path:
        p_proj, proj_mu = ObservableProj(p,mu_old)
        proj.append(p_proj)
        mu_new.update(proj_mu)
    T_new = proj

    #marcando estados
    Xm_new = []
    for x in X_new:
        if x in G.Xm:
            Xm_new +=[x]
    
    # Criando o autômato projeção e juntando com a função de tempo mu_new
    Pi_Gt = fsa(X_new, E_new, T_new, X0_new, Xm_new, name = '$\Pi(' + G.name + ')$')
    Pi_Gt.setgraphic(style='observer')
    Gproj = (Pi_Gt,mu_new)
                    
    return Gproj


## Retorna a interseção entre todos os intervalos.
def mu_intersection(mu, intervals):
    inter = P.closedopen(0,P.inf)
##    print("o que é isso: " + str(intervals))
    for i in intervals:
##        print("isso que esta tentando fazer a operação: " + str(i))
        inter = inter & i
        
    return inter

## Constrói o automato determinístico equivalente
def ti_equi_det(Gt):
    G = Gt[0]
    mu_old = Gt[1]
    initial = list(G.X0)
    if len(initial) > 1:
        X0_new = tuple(initial)
        X_new = [deepcopy(X0_new)]
    else:
        X0_new = list(G.X0)
        X_new = X0_new.copy()
    Xm_new = [] 
   
    #X_new = [X0_new]
    E_new = list(G.Sigma)
    T_new = []
    mu_new = dict()

    X_aux = X_new.copy()
    
    while len(X_aux) != 0:
        x_actual = X_aux.pop(0)
        
        if  isinstance(x_actual, tuple):
            ## Recuperar transições de todos os estados de dentro
            
            trans_dict = {t:i for (t,i) in mu_old.items() if t[0] in x_actual}
                
            active_events = set()
            
            for q in x_actual:
                active_events = active_events | G.Gamma(q) #uniao de conjuntos
            #active_events = active_events & set(E_new) #intersecao de conjuntos

            for ev in active_events:
                t_aux = [t for (t,i) in trans_dict.items() if t[1] == ev]
                nm,d = max_disj_trans(mu_old,t_aux) #conjuntos disjuntos

                nm_upd = {(x_actual,sig,x_dest):i for (x_orig,sig,x_dest),i in nm.items()}
                mu_new.update(nm_upd)

                T_new += [t for (t,i) in nm_upd.items()]
                
                x_dummy = [t[2] for (t,i) in nm_upd.items() if t[2] not in X_new]
                x_dummy = sorted(x_dummy, key=Py2Key)
                X_new += x_dummy
                X_aux += x_dummy
            
        else: 

            trans_dict = dict()

            #pega todas as trasicoes que saem do estado no automato original
        
            trans_dict = {t:i for (t,i) in mu_old.items() if x_actual == t[0]}
       
            for ev in G.Gamma(x_actual):
                
                t_aux = [t for (t,i) in trans_dict.items() if t[1] == ev]
                nm,d = max_disj_trans(mu_old,t_aux)

                mu_new.update(nm)
                T_new += [t for (t,i) in nm.items()]
                
                x_dummy = [t[2] for (t,i) in nm.items() if t[2] not in X_new]
                x_dummy = sorted(x_dummy, key=Py2Key)
                X_new += x_dummy
                X_aux += x_dummy


    #definindo estados marcados
    for x in X_new:
        if isinstance(x,tuple):
            xlist = list(x)
            if len(list(set(xlist)& G.Xm)) != 0:
                Xm_new += [x]
        else:
            if x in G.Xm:
                Xm_new += [x]

                
    Gdet = fsa(X_new,E_new,T_new,X0_new,Xm_new,name = "$G_{det}$")
    Gdet.setgraphic(style='observer')
    Gt_det = tia(Gdet,mu_new)
    
    return Gt_det


## Retorna um dicionario com as transições disjuntas, incluindo 
## novos estados e seus I associados.
def max_disj_trans(mu,transitions):

    ## ERROR HANDLE ZONE ===============================================\/
##    print("O que esta aparecendo em *transitions? R: " + str(transitions))
    result = all(sig == transitions[0][1] for (x,sig,y) in transitions)
    if not result:
##        print("Todas as transições são rotuladas pelo mesmo evento.")
##    else:
        raise NameError("All transitions must be labeled by the same event!")

    ## =================================================================/\
    from itertools import combinations as comb
    
    disj_set = []

    oldmu = {t:i for (t,i) in mu.items() if t in transitions}
    newmu = dict()
    
    t_int = [oldmu[t] for t in transitions]
##    print("todos os intervalos disponíveis: " + str(t_int))
##    print("--------------------------------")
##    print("")
    
##    for i in range(len(t_int)):
    for i in range(len(t_int),0,-1):
##        print("--------------------------")
##        print("combinacoes de " + str(i) + " intervalos")
        aux = []
        for c in comb(t_int,i):
            aux.extend(mu_intersection(mu,c))
##        print("")
##        print("intersecoes: " + str(aux))
##        print("")
        for inter in range(len(aux)):
            for disj_int in disj_set:
##                print(str(aux[inter]) + " - " + str(disj_int))
                aux[inter] = aux[inter] - disj_int
##                print("= " + str(aux[inter]))
##                if aux[inter] == P.empty():
##                    aux.pop(inter)
        disj_set.extend([it for it in aux if it != P.empty()])
##        print("intervalos disj calculados: " + str(aux))
##        disj_set.extend(aux)
        
##    print("==========================")
    disj_set.sort()
##    print("intervalos disj calculados: " + str(disj_set))
##    print("--------------------------------")
##    print("")

    # TODO: 
    #       No fim da comparação, cria a transição do estado original e evento
    #       para o estado que contenha todos os estados gravados.
    #       Se a transição ja existir, faz união com o intervalo ja existente.
    
    for disj in disj_set:
        x_aux = ()
        for orig in oldmu:
            if disj & oldmu[orig]:
##                print(str(disj) + " faz interseção com " + str(oldmu[orig]) + ", o estado destino é: " + str(orig[2]))
                x_aux = x_aux + tuple([orig[2]])
##            else:
##                print(str(disj) + " NÃO tem interseção com " + str(oldmu[orig]) + "!")
        if len(x_aux) > 1:
            xa = list(x_aux)
            xa = sorted(xa, key=Py2Key)
##            xa.sort()
            
            trans_aux = (orig[0],orig[1],tuple(xa))
        else:
            trans_aux = (orig[0],orig[1],x_aux[0])
##        print("")
##        print("Transição para o intervalo " + str(disj) + " é " + str(trans_aux))
##        print("")
        if trans_aux in newmu:
            newmu.update({trans_aux : newmu[trans_aux]|disj})
        else:
            newmu.update({trans_aux : disj})
##        print("Novo mu: " + str(newmu))
##        print("--------------------------------")
##        print("")
    
    return newmu, disj_set


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
                                    tr_from = (xr,t[1],t[2])
                                    Tr.remove(t)

                                    if p[-1][2] in Xr:
##                                        print("Removendo estado " + str(p[-1][2]))
                                        Xr.remove(p[-1][2])
                                    Tr.append(tr_from)

                                    Murf[tr_from] = Mu[t]
                                    Murf.pop(t)
##                                    print("Renamed Transition From: " + str(tr_from))
                                

    ##                    tr_from = [((xo,evr,Mu[(xo,e,xd)]),e,xd) for (xo,e,xd) in  if xo==p[-1][2]]

                        Xr.append(xr)

                    elif ret == 'GRI':
                        tr_to = (p[-1][0],evr,p[-1][2])

##                    print("Removendo transição " + str(p[-1]))
                    Tr.remove(p[-1])
                    Tr.append(tr_to)

##                    for t in tr_from:
##                        Tr.remove
                    if ret == 'GRI':
                        Muri[tr_to] = Mu[p[-1]]
                        Muri.pop(p[-1])
                    elif ret == 'GRF':
                        Murf[tr_to] = Mu[p[-1]]
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
        glrf = ac(fsa(Xr,Er,Tr,X0r,Sigobs = Eor,name="$G_{lT}^R$"))
        glrf.setgraphic(style='observer')
        glt_rf = tia(glrf,Murf)
        return glt_rf

## Em Construção: Descompacta os estados renomeados de Gdtr,
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
            
            xunp.sort()
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

## Realiza a operação produto entre dois TIA
def ti_product(self, other):
    
    def pars(var):
        resp=()
        for i in var:
            if type(i)==tuple:
                resp+=pars(i)
            else:
                resp+=(i,)
        return resp
    
    def Gamma_prod(p):
        """ definition of the Gamma function of parallel
        composition of the input automata for the tuple
        p = (x1,x2)"""
        x1, x2 = p
        active_events = Gamma1(x1) & Gamma2(x2)
        
        return active_events

    def latexname(p,simplify):  
        """This function takes a tuple of the composite 
        automaton and generate a latex pretty name, based
        on the dictionaries of each input automaton"""
        x1, x2 = p
        if x1 in self[0].symDict:
            x1_name = self.symDict[x1]
        else: 
            x1_name = str(x1)
        if x2 in other[0].symDict:
            x2_name = other.symDict[x2]
        else: 
            x2_name = str(x2)
        name = '('+x1_name+','+x2_name+')'
        if simplify:
            name = name.replace("(","")
            name = name.replace(")","")
            name = '('+name+')'
        return name

    # Main program

    # G1 marks an empty language and G2 does not
    if self[0].empty and not other[0].empty:
        return tia(fsa(),dict())
    # G1 and G2 marks empty languages
    elif self[0].empty and other[0].empty:
        return tia(fsa(),dict())
    # G1 marks a nonempty language and G2 does.
    elif not self[0].empty and other[0].empty:
        return tia(fsa(),dict())
    
    # easy names for inner variables
    X1, X2 = self[0].X, other[0].X
    X01, X02 = self[0].X0, other[0].X0
    Xm1, Xm2 = self[0].Xm, other[0].Xm 
    Sigma1, Sigma2 = self[0].Sigma, other[0].Sigma
    Gamma1, Gamma2 = self[0].Gamma, other[0].Gamma 
    delta1, delta2 = self[0].delta, other[0].delta   
    Sigobs1, Sigobs2 = self[0].Sigobs, other[0].Sigobs         
    Sigcon1, Sigcon2 = self[0].Sigcon, other[0].Sigcon
    Mu1, Mu2 = self[1], other[1]

    # initial values for composite automaton
    X0_p = [(x01,x02) for x01 in X01 for x02 in X02]
    S = [(x01,x02) for x01 in X01 for x02 in X02]
    X_p = set(X0_p)
    Xm_p = frozenset([])
    Sigobs_p = Sigobs1|Sigobs2 
    Sigcon_p = Sigcon1|Sigcon2
    Sigma_p = Sigma1 | Sigma2
    transition=[]
    mu_p = dict()

    
    # generating labeling table for initial states
    table_p = {}
    for p in X0_p: # labeling initial states
        #table_p.update({p: latexname(p,simplify)})
        p1,p2 = p
        if (p1 in Xm1) & (p2 in Xm2):
            Xm_p = Xm_p | set([p])
   
    # main stack cycle    
    while S != [] :
        p = S.pop() # taking a tuple
        #lembrando que p1 e p2 podem ser tuplas
        p1, p2 = p     
        for sigma in Gamma_prod(p):
            Q1_n, Q2_n = delta1(p1, sigma), delta2(p2,sigma)
            if type(Q1_n) == tuple:
                Q1_n = set([Q1_n])
            if type(Q2_n) == tuple:
                Q2_n = set([Q2_n])
##            print("Q1n: " + str(Q1_n) + " " + str(type(Q1_n)) + "\nQ2n: " + str(Q2_n) +" " + str(type(Q2_n)) +"\n")
            try:
                type1 = (type(Q1_n) == str or type(Q1_n) == int)
                type2 = (type(Q2_n) == str or type(Q2_n) == int)
                if type1 and not type2:
##                    print("Q1 é str e Q2 não!")
                    Q = list()
                    for q2 in Q2_n:
                        Q.append((Q1_n,q2))
                elif not type1 and type2:
##                    print("Q2 é str e Q1 não!")
                    Q = list()
                    for q1 in Q1_n:
                        Q.append((q1,Q2_n))
                elif type1 and type2:
##                    print("Tanto Q1 quanto Q2 são str!")
                    #Q = list((Q1_n,Q2_n))
                    Q = [(Q1_n,Q2_n)]
                elif not type1 and not type2:
##                    print("Nem Q1 nem Q2 são str!")
                    #TODAS AS DUPLAS POSSIVEIS
                    Q = list((q1,q2) for q1 in Q1_n  for q2 in Q2_n)
##                print("Q: " + str(Q) + "\n")
                
                
                # Q contains the tuples for reached states
                for q in Q:
##                    print("p1: " + str(p1) + " | p2: " + str(p2))
##                    print("q: " + str(q) + "\n")
                    if Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])] != P.empty():
                        if q not in X_p:
                            # updating states of composite automaton and stack
                            X_p = X_p | set([q])           
                            S.append(q)
                            q1, q2 = q
                            # updating marked states
                            if (q1 in Xm1) & (q2 in Xm2):
                                Xm_p = Xm_p | set([q]) 
                            #updating latex dictionary
                            #table_p.update({q: latexname(q,simplify)})
                        transition.append([p,sigma,q])
                        if not mu_p.get((p,sigma,q)):
                            mu_int = Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])]
                            mu_p[(p,sigma,q)] = mu_int
##                print("--------------\n")
            except:
##                print("> EXCEPT <")
                q = (Q1_n,Q2_n)
##                print("p1: " + str(p1) + " | p2: " + str(p2))
##                print("q: " + str(q) + "\n")
##                print("--------------\n")
                if Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])] != P.empty():
                    if q not in X_p:
                        # updating state and stack                                                     
                        X_p = X_p | set([q])           
                        S.append(q)
                        q1, q2 = q
                        if (q1 in Xm1) & (q2 in Xm2):
                            # updating marked states
                            Xm_p =Xm_p | set([q])
                        # updating labeling table 
                        #table_p.update({q: latexname(q,simplify)})
                    transition.append([p,sigma,q])
                    if not mu_p.get((p,sigma,q)):
                        mu_int = Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])]
                        mu_p[(p,sigma,q)] = mu_int
                    
            #if sigma in self[0].symDict:
                #table_p.update({sigma: self[0].symDict[sigma]}) 
            #elif sigma in other[0].symDict:
                #table_p.update({sigma: other[0].symDict[sigma]}) 

    # finally, we build the composite automaton    
   # prodnondet = fsa(list(X_p), Sigma_p, transition, X0_p, list(Xm_p), table=table_p,
                  #Sigobs = Sigobs_p, Sigcon=Sigcon_p,name='untitled')
    prodnondet = fsa(list(X_p), Sigma_p, transition, X0_p, list(Xm_p), table=table_p,
                  Sigobs = Sigobs_p, Sigcon=Sigcon_p,name='untitled')
    tiprod = tia(prodnondet,mu_p)
    return tiprod

## Faz o produto extendido entre dois TIA, colocando self loop com
## eventos particulares do outro autômato em todos os estados
def ext_ti_product(self,other):
    G1, Mu1 = self
    G2, Mu2 = other    
    sigp_1 = G1.Sigma - G2.Sigma
    sigp_2 = G2.Sigma - G1.Sigma
##    print("Sigma particular de G1: " + str(sigp_1))
##    print("Sigma particular de G2: " + str(sigp_2))
    
    for sig in sigp_2:
        for x in G1.X:
            G1 = G1.addselfloop(x,sig)
            Mu1[(x,sig,x)] = P.closedopen(0,P.inf)
##            print("Adicionando self loop de " + str(sig) + " no estado " + str(x) + " de G1")
##    print("\n Sigma de G1 agora: " + str(G1.Sigma))
    for sig in sigp_1:
        for x in G2.X:
            G2 = G2.addselfloop(x,sig)
            Mu2[(x,sig,x)] = P.closedopen(0,P.inf)
##            print("Adicionando self loop de " + str(sig) + " no estado " + str(x) + " de G2")
##    print("\n Sigma de G2 agora: " + str(G2.Sigma))

    extprod = ti_product(tia(G1,Mu1),tia(G2,Mu2))
    extprod[0].setgraphic(style='observer')
    return extprod
        
    
## Adiciona os tempos maximo e mínimo de dois intervalos
def interval_offset(self,other):    
    left1, left2 = self.left , other.left
    lower1, lower2 = self.lower , other.lower
    upper1, upper2 = self.upper , other.upper
    right1, right2 = self.right , other.right
    if left1 == P.OPEN or left2 == P.OPEN:
        left_off = False
    else:
        left_off = True
    lower_off = lower1 + lower2
    upper_off = upper1 + upper2
    if right1 == P.OPEN or right2 == P.OPEN:
        right_off = False
    else:
        right_off = True
    return P.from_data([(left_off,lower_off,upper_off,right_off)])

## Em Construção: mesma função simplify da toolbox Diagnosis, mas para um TIA
def ti_simplify(Gt):
    g=Gt[0].copy()
    mu=Gt[1].copy()
    renmu = dict()
    mapping=[]
    for state in g.X:
        tex="".join(str(i) for i in state)
        tex=tex.replace('(',"")
        tex=tex.replace(')',"")
        tex=tex.replace(',',"")
        tex=tex.replace(' ',"")
        tex=tex.replace("'","")
        mapping.append((state,tex))
    for item in mu.items():
        t,interval = item
        x,sig,y = t
        newx = str(x)
        newx = newx.replace('(',"")
        newx = newx.replace(')',"")
        newx = newx.replace(',',"")
        newx = newx.replace(' ',"")
        newx = newx.replace("'","")
        newy = str(y)
        newy = newy.replace('(',"")
        newy = newy.replace(')',"")
        newy = newy.replace(',',"")
        newy = newy.replace(' ',"")
        newy = newy.replace("'","")
        renmu[(newx,sig,newy)]=interval
##        print(renmu)
    renamed_gt = tia(g.renamestates(mapping),renmu)
    return renamed_gt

## Realiza o processo de construção dos automatos para o diagnóstico
## de falhas em TIA.
##ret = "GD" : (padrão) Retorna o Diagnosticador usando ti proj
##ret = "AL" : Retorna o automato rotulador com intervalo de tempo
##ret = "GL" : Retorna o produto de Gt com Alt, automato rotulado
def ti_diag(Gt,failevent,ret="GD"):
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
    
    Glt = ti_product(Gt,Alt)
    Glt[0].setgraphic(style='observer')

    if ret == "GL":
        return Glt
    else:
        pi_glt = ti_proj(Glt)
        diag = ti_equi_det(pi_glt)
        return diag

## Retorna True/False caso o Gt tenha SCC formados por estados incertos
## TODO: NÃO ESTÁ FUNCIONANDO TOTALMENTE... precisa ajeitar para todos os casos
def has_scc(Gt, ret = False):
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


#Complementar de um TIA
def ti_complement(Gt):
    #desmarcar todos os estados e adicionar estado dump xd
    G = Gt[0]
    mu_new = Gt[1].copy()
    NewMarked = list(G.X - G.Xm) + ['dump']
    #NewStates = list(G.X) + ['dump']
    G = G.addstate('dump')
    G = G.setpar(Xm=NewMarked)
    
    Gtcomp = tia(G,mu_new)

    for x in G.X:
        trans = []
        active = []
        #conjunto de transicoes que saem de x
        for t in G.transitions():
            if t[0]==x:
                trans += [t]
                active += [t[1]] #eventos ativos
        #fazer o complementar da uniao dos intervalos associados ao evento ev
        #que sai do estado x
        for ev in active:
            union = P.empty()
            for t in trans:
                if t[1] == ev:
                    union = union | mu_new[t]
            comp = P.closedopen(0,P.inf) - union #intervalo complementar
            #criar transicao com evento ev e intervalo comp que leva para dump
            newt = (x,ev,'dump') #nova transicao
            mu_new[newt] = comp #intervalo associado a nova transicao
            G = G.addtransition(newt)
            Gtcomp = tia(G, mu_new)

        #transicoes para o dump com os eventos nao definidos no estado 
        for ev in list(set(G.Sigma) - set(active)):
            newt = (x,ev,'dump')
            mu_new[newt] = P.closedopen(0,P.inf)
            G = G.addtransition(newt)
            Gtcomp = tia(G,mu_new)

    #autolacos no dump
    for ev in list(G.Sigma):
        newt = ('dump',ev,'dump')
        mu_new[newt] = P.closedopen(0,P.inf)
        G = G.addtransition(newt)
        Gtcomp = tia(G,mu_new)
        

    return Gtcomp
            

#-----------------------------------------------------------------------------
#-------------------------- opacidade novo -----------------------------------
#-----------------------------------------------------------------------------

## Realiza a operação produto entre dois TIA colocando labels
def ti_label_obf(self, other):
    
    def pars(var):
        resp=()
        for i in var:
            if type(i)==tuple:
                resp+=pars(i)
            else:
                resp+=(i,)
        return resp
    
    def Gamma_prod(p):
        """ definition of the Gamma function of parallel
        composition of the input automata for the tuple
        p = (x1,x2)"""
        x1, x2 = p
        active_events = Gamma1(x1) & Gamma2(x2)
        
        return active_events

    def latexname(p,simplify):  
        """This function takes a tuple of the composite 
        automaton and generate a latex pretty name, based
        on the dictionaries of each input automaton"""
        x1, x2 = p
        if x1 in self[0].symDict:
            x1_name = self.symDict[x1]
        else: 
            x1_name = str(x1)
        if x2 in other[0].symDict:
            x2_name = other.symDict[x2]
        else: 
            x2_name = str(x2)
        name = '('+x1_name+','+x2_name+')'
        if simplify:
            name = name.replace("(","")
            name = name.replace(")","")
            name = '('+name+')'
        return name

    # Main program

    # G1 marks an empty language and G2 does not
    if self[0].empty and not other[0].empty:
        return tia(fsa(),dict())
    # G1 and G2 marks empty languages
    elif self[0].empty and other[0].empty:
        return tia(fsa(),dict())
    # G1 marks a nonempty language and G2 does.
    elif not self[0].empty and other[0].empty:
        return tia(fsa(),dict())
    
    # easy names for inner variables
    X1, X2 = self[0].X, other[0].X
    X01, X02 = self[0].X0, other[0].X0
    Xm1, Xm2 = self[0].Xm, other[0].Xm 
    Sigma1, Sigma2 = self[0].Sigma, other[0].Sigma
    Gamma1, Gamma2 = self[0].Gamma, other[0].Gamma 
    delta1, delta2 = self[0].delta, other[0].delta   
    Sigobs1, Sigobs2 = self[0].Sigobs, other[0].Sigobs         
    Sigcon1, Sigcon2 = self[0].Sigcon, other[0].Sigcon
    Mu1, Mu2 = self[1], other[1]


    # initial values for composite automaton
    X0_p = [(x01,x02) for x01 in X01 for x02 in X02]
    S = [(x01,x02) for x01 in X01 for x02 in X02]
    X_p = set(X0_p)
    Xm_p = frozenset([]) #DEPOIS TRANSFORMAR EM LISTA?
    #Sigobs_p = Sigobs1|Sigobs2 
    #Sigcon_p = Sigcon1|Sigcon2
    Sigma_p = list(Sigma1 | Sigma2)
    Sigma_pl = []
    transition=[]
    mu_p = dict()

    #renaming events
    for el in Sigma_p:
        Sigma_pl.append(el+' '+'co')
        Sigma_pl.append(el+' '+'po')
        Sigma_pl.append(el+' '+'cr')
        Sigma_pl.append(el+' '+'pr')

    # generating labeling table for initial states
    table_p = {}
    for p in X0_p: # labeling initial states
        #table_p.update({p: latexname(p,simplify)})
        p1,p2 = p
        if (p1 in Xm1) & (p2 in Xm2):
            Xm_p = Xm_p | set([p])
   
    # main stack cycle    
    while S != [] :
        p = S.pop() # taking a tuple
        #lembrando que p1 e p2 podem ser tuplas
        p1, p2 = p     
        for sigma in Gamma_prod(p):
            Q1_n, Q2_n = delta1(p1, sigma), delta2(p2,sigma)
            if type(Q1_n) == tuple:
                Q1_n = set([Q1_n])
            if type(Q2_n) == tuple:
                Q2_n = set([Q2_n])
##            print("Q1n: " + str(Q1_n) + " " + str(type(Q1_n)) + "\nQ2n: " + str(Q2_n) +" " + str(type(Q2_n)) +"\n")
            try:
                type1 = (type(Q1_n) == str or type(Q1_n) == int)
                type2 = (type(Q2_n) == str or type(Q2_n) == int)
                if type1 and not type2:
##                    print("Q1 é str e Q2 não!")
                    Q = list()
                    for q2 in Q2_n:
                        Q.append((Q1_n,q2))
                elif not type1 and type2:
##                    print("Q2 é str e Q1 não!")
                    Q = list()
                    for q1 in Q1_n:
                        Q.append((q1,Q2_n))
                elif type1 and type2:
##                    print("Tanto Q1 quanto Q2 são str!")
                    #Q = list((Q1_n,Q2_n))
                    Q = [(Q1_n,Q2_n)]
                elif not type1 and not type2:
##                    print("Nem Q1 nem Q2 são str!")
                    #TODAS AS DUPLAS POSSIVEIS
                    Q = list((q1,q2) for q1 in Q1_n  for q2 in Q2_n)
##                print("Q: " + str(Q) + "\n")
                
                
                # Q contains the tuples for reached states
                for q in Q:
##                    print("p1: " + str(p1) + " | p2: " + str(p2))
##                    print("q: " + str(q) + "\n")
                    if Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])] != P.empty():
                        if q not in X_p:
                            # updating states of composite automaton and stack
                            X_p = X_p | set([q])           
                            S.append(q)
                            q1, q2 = q
                            # updating marked states
                            if (q1 in Xm1) & (q2 in Xm2):
                                Xm_p = Xm_p | set([q]) 
                            #updating latex dictionary
                            #table_p.update({q: latexname(q,simplify)})
                                
                        if Mu1[(p1,sigma,q[0])] in Mu2[(p2,sigma,q[1])]:
                            transition.append([p,sigma+' '+'co',q])
                            if not mu_p.get((p,sigma+' '+'co',q)):
                                mu_int = Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])]
                                mu_p[(p,sigma+' '+'co',q)] = mu_int    
                        else:
                            transition.append([p,sigma+' '+'po',q])
                            if not mu_p.get((p,sigma+' '+'po',q)):
                                mu_int = Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])]
                                mu_p[(p,sigma+' '+'po',q)] = mu_int    

                           
##                print("--------------\n")
            except:
##                print("> EXCEPT <")
                q = (Q1_n,Q2_n)
##                print("p1: " + str(p1) + " | p2: " + str(p2))
##                print("q: " + str(q) + "\n")
##                print("--------------\n")
                if Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])] != P.empty():
                    if q not in X_p:
                        # updating state and stack                                                     
                        X_p = X_p | set([q])           
                        S.append(q)
                        q1, q2 = q
                        if (q1 in Xm1) & (q2 in Xm2):
                            # updating marked states
                            Xm_p =Xm_p | set([q])
                        # updating labeling table 
                        #table_p.update({q: latexname(q,simplify)})

                    if Mu1[(p1,sigma,q[0])] in Mu2[(p2,sigma,q[1])]:
                        transition.append([p,sigma+' '+'co',q])
                        if not mu_p.get((p,sigma+' '+'co',q)):
                            mu_int = Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])]
                            mu_p[[p,sigma+' '+'co',q]] = mu_int    
                        else:
                            transition.append([p,sigma+' '+'po',q])
                            if not mu_p.get((p,sigma+' '+'po',q)):
                                mu_int = Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])]
                                mu_p[[p,sigma+' '+'po',q]] = mu_int   

                     #transition.append([p,sigma,q])

                   
                        #mu_p[(p,sigma,q)] = mu_int
                    
            #if sigma in self[0].symDict:
                #table_p.update({sigma: self[0].symDict[sigma]}) 
            #elif sigma in other[0].symDict:
                #table_p.update({sigma: other[0].symDict[sigma]}) 

    # finally, we build the composite automaton    
   # prodnondet = fsa(list(X_p), Sigma_p, transition, X0_p, list(Xm_p), table=table_p,
                  #Sigobs = Sigobs_p, Sigcon=Sigcon_p,name='untitled')
    prodnondet = fsa(list(X_p), Sigma_pl, transition, X0_p, list(Xm_p), table=table_p,
                     name='$G_{obf}$')
    tiprod = tia(prodnondet,mu_p)
    return tiprod

    
    


## Realiza a operação produto entre dois TIA colocando labels Grev
#Gs x Gns^c

def ti_label_rev(self, other):
    othercomp = ti_complement(other)
    Grev = ti_label_obf(self,othercomp)
    Comp_marked = othercomp[0].Xm
    G = Grev[0]
    dic = Grev[1] 
    #M = Grev[0].X

    for t in G.transitions():
        x = t[0]
        ev = t[1]
        label = ev[len(ev)-2:] #label po, co pega os dois ultimos
        event = ev[:-3] #evento em si, a, b, c ... tira os dois ultimos
        xprox = t[2]
        time = dic[t]
        #if xprox in M:
        if (x[1] != 'dump' and (xprox[1] in Comp_marked)):
            if label == 'po':
                new = event+' '+'pr'
            if label == 'co':
                new = event+' '+'cr'
                    
            G = G.renametransition([x,(ev,new),xprox])
            dic[(x,new,xprox)] = time

    Grev = fsa(G.X, G.Sigma, G.transitions(), G.X0, G.Xm,
                     name='$G_{rev}$')
    Grev = tia(Grev,dic)
    return Grev
            
        
def verifierTLBO(secret,nsecret,self,other):
    #automatos
    #self = Gobf
    #other = Grev
    #secret = Gs
    #nsecret = Gns

    #dicionarios:
    #dself = label_obf
    #dother = label_rev

    dself = {} #dicionario com os labels {transicao - label}
    dother = {}
    self = self[0]
    other = other[0]

    for t in self.transitions():
        ev = t[1]
        e = ev[:-3] #tem espaço, por exemplo 'a co' entao tem que tirar os 3 ultimos pra pegar o nome do evento
        lbl = ev[-2:]
        self = self.renametransition([t[0],(ev,e),t[2]]) #tirando label do evento
        #Gobf[1][t[0],e,t[2]] = Gobf[1].pop[t[0],ev,t[2]] nao preciso mais dos tempos entao talvez essa linha seja inutil
            
        dself[(t[0],e,t[2])] = lbl #adicionando label ao dicionario

    for t in other.transitions():
        ev = t[1]
        e = ev[:-3]
        lbl = ev[-2:]
        other = other.renametransition([t[0],(ev,e),t[2]]) #tirando label do evento

        dother[(t[0],e,t[2])] = lbl
        

    newdict = {} #novo dicionario para os labels das transicoes
    #automato verifier vai ser o fsa + newdict = "tia" (acho que da pra fazer assim) 

        
    # -----------------------------------------------------------------------------------------------------------
        
    # Main program 
    # easy names for inner variables
    X1, X2 = self.X, other.X
    X01, X02 = self.X0, other.X0
    Xm1, Xm2 = self.Xm, other.Xm 
    Sigma1, Sigma2 = secret.Sigma, nsecret.Sigma 
    Gamma1, Gamma2 = self.Gamma, other.Gamma 
    delta1, delta2 = self.__delta__, other.__delta__    
    #Sigobs1,Sigobs2 = secret.Sigobs, nsecret.Sigobs         
    #Sigcon1,Sigcon2 = secret.Sigcon, nsecret.Sigcon  
        
    # initial values for the composite automaton
    #X0_p = [(tuple(X01),tuple(X02))] #Gobf e Grev normalmente nem tem mais de um estado inicial
    X0_p = []
    state_dict = {}

    for state in secret.X0:
        init_obf = ()
        init_rev = ()
        for es_obf in X01:
            if es_obf[0] == state:
                init_obf += (es_obf,)
        for es_rev in X02:
            if es_rev[0] == state:
                init_rev += (es_rev,)
        state_dict[(init_obf,init_rev)] = state
        X0_p += [(init_obf,init_rev)]
        
            
    S = deepcopy(X0_p)
    X_p = deepcopy(X0_p)
    Xm_p = []
    #Sigobs_p = Sigobs1|Sigobs2 
    #Sigcon_p = Sigcon1|Sigcon2
    Sigma_p = Sigma1 | Sigma2
    transition=[]

    # main stack cycle     
    while S != [] :
        p = S.pop(); # taking a tuple
        state_p = state_dict[p]
        #print(state_p)
        es1 = p[0]
        es2 = p[1]
        #es1 e es2 podem ser conjuntos de tuplas ou so uma tupla        

        #transformei tudo em lista de uma vez
        #acho que nao precisa, tudo parece ja estar setado como lista mesmo
        if type(es1) != tuple:
            es1 = tuple(es1)
        if type(es2) != tuple:
            es1 = tuple(es2)

        for ev in Sigma1 | Sigma2: #aqui talvez seja melhor só escanear pelo conjunto de eventos ativos em ambos
            Next1 = []
            Next2 = []

            #listas de proximos estados executando o evento ev
            for state1 in es1:
                if state1 != 'empty':
                    if ev in Gamma1(state1):
                        Next1 += tuple(delta1(state1,ev))
                if state1 == 'empty':
                    Next1 = ('empty',) #se for empty continua, nao soma porque quero um empty só
            for state2 in es2:
                if state2 != 'empty':
                    if ev in Gamma2(state2):
                        Next2 += tuple(delta2(state2,ev))
                if state2 == 'empty':
                    Next2 = ('empty',) #se for empty continua, nao soma porque quero um empty só


            
                
            #checar estados de Gns na lista de proximos estados de Gobf e Grev
            #nao sei se essa é a forma mais eficiente ... 
            for state in secret.X:
                lab_obf = []
                lab_rev = []
                
                Next1_state = () #lista de proximos estados de Gobf cuja primeira componente é o estado state em Gs
                Next2_state = () #lista de proximos estados de Grev cuja primeira componente é o estado state em Gs

                #se um deles ja for empty, aqui vai virar a tupla vazia ()
                for es in Next1:
                    if es[0] == state:
                        Next1_state += (es,)
                    if type(es[0]) == tuple:
                        if state in es[0]:
                            Next1_state += (es,)
                for es in Next2:
                    if es[0] == state:
                        Next2_state += (es,)
                    if type(es[0]) == tuple:
                        if state in es[0]:
                            Next2_state += (es,)


                #os dois nao podem ser vazios ao mesmo tempo. se forem, o automato nao anda pra frente
                #aqui funciona tanto para o caso de nao ter proxima transicao quanto para o caso de ja ser empty antes (propagar o empty)
                #só transformo em empty quando um deles nao é, entao se os dois forem vazios vai continuar [] para ambos
                if Next1_state == () and Next2_state != ():
                    Next1_state = ('empty',)
                if Next2_state == () and Next1_state != ():
                    Next2_state = ('empty',)

                if Next1_state != () or Next2_state != (): #os dois nao podem ser vazios ao mesmo tempo
                    #new_state_double = (Next1_state, Next2_state) #novo estado, pode ter duplicatas em cada componente
                    
                    #removendo redundancias (estados duplicados)
                    Next1_state = tuple(set(Next1_state))
                    Next2_state = tuple(set(Next2_state))

                    new_state = (Next1_state, Next2_state)
                    state_dict[new_state] = state
                    #print(ev)
                    #print(state)

                    if (state_p,ev,state) in secret.transitions():
                        new_transition = (p,ev,new_state)
              
                        if new_transition not in transition:
                            #adicionar estado e transicao ao automato e a lista S tambem
                            S.append(new_state)
                            X_p.append(new_state)
                            transition.append(new_transition)

                            #estados marcados
                            if state in secret.Xm:
                                Xm_p.append(new_state)

                            #adicionar tudo nos dicionarios e no dicionario novo
                            for pobf in p[0]:
                                for next_obf in Next1_state:
                                    tobf = (pobf,ev,next_obf)
                                    if tobf in self.transitions():
                                        lab_obf += [dself[tobf]]

                            for prev in p[1]:
                                for next_rev in Next2_state:
                                    trev = (prev,ev,next_rev)
                                    if trev in other.transitions():
                                        lab_rev += [dother[trev]]

                            if lab_obf == []:
                                lab_obf = ['empty']
                            if lab_rev == []:
                                lab_rev = ['empty']

                            #acho que vou ter que trocar tudo pra tupla (pelo menos os objetos, que sao new_transition)
                            newdict[new_transition] = [lab_obf,lab_rev] #dicionario novo
                       

                        
    # finally, we build the composite automaton                                
    verifier = fsa(X_p, Sigma_p, transition, X0_p, Xm_p,name='Verifier')
    verifier.setgraphic(style='rectangle')
    #verifier_label = tia(verifier,newdict)

    return verifier, newdict, state_dict


def label_draw(*Gt):
    if type(Gt[-1]) == str:
        style = Gt[-1]
    else:
        style = 'beamer'
    for tia in Gt:
        if type(tia) != tuple:
            break
        Gtia = tia[0]
        new_sigobs = []
        new_sig = []
        if not isitempty(Gtia):
            for trans in transitions(tia[0]):
                ev_mu = trans[1]+str(tia[1][trans])
                if trans[1] in tia[0].Sigobs:
                    new_sigobs.append(ev_mu)
                new_sig.append(ev_mu)
                Gtia = Gtia.renametransition([trans[0],(trans[1],ev_mu),trans[2]])
            Gtia = Gtia.setpar(Sigobs = new_sigobs)
            draw(Gtia,style)
        else:
            draw(fsa())
    return






    
