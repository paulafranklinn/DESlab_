# Functions for Time-Interval Discrete Event Systems

from deslab import *
import portion as P
from sys import intern
from copy import deepcopy

## Documents for this library can be found here:
## https://github.com/AlexandreDecan/portion

#AUXILIARY CLASSES:
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
    """
    This instruction defines a time-interval automaton object inside DESlab.

    -------
    Example

    syms('a b u')
    # automaton definition G_T
    Xt = [0, 1, 2, 3, 4, 5, 6]
    Et =[a,b,u]
    sigobst = [a,b,u]
    X0t = [0]
    Xmt = [ ]
    Tt = [(0,u,1),(0,a,4),(1,a,2),(2,b,3),(4,u,5),(5,b,6),(5,b,3)]
    
    mut = {(0,u,1): P.closed(0,2),
          (0,a,4): P.closed(1,3.5),
          (1,a,2): P.closed(1,3),
          (2,b,3): P.closed(2,3),
          (4,u,5): P.closed(0.5,1.5),
          (5,b,6): P.closed(2,4),
          (5,b,3): P.closed(1,3)
          }
    
    G = fsa(Xt,Et,Tt,X0t,Xmt,Sigobs=sigobst,name="$G_T$")
    GT = tia(G,mut)
    ti_draw(GT)
  
    """
    Gt = (G,mu)
    return Gt


def ti_draw(*Gt):
    """
    This instruction Create the state transition diagram of a 
    time-interval automaton where transitions consist of the event 
    and the corresponding time interval.
    It receives the TIA and if desired the display mode:
        Figure
        Figurecolor
        Beamer
    
    -------
    Example

    syms('a b u')
    # automaton definition G_T
    Xt = [0, 1, 2, 3, 4, 5, 6]
    Et =[a,b,u]
    sigobst = [a,b,u]
    X0t = [0]
    Xmt = [ ]
    Tt = [(0,u,1),(0,a,4),(1,a,2),(2,b,3),(4,u,5),(5,b,6),(5,b,3)]
    
    mut = {(0,u,1): P.closed(0,2),
          (0,a,4): P.closed(1,3.5),
          (1,a,2): P.closed(1,3),
          (2,b,3): P.closed(2,3),
          (4,u,5): P.closed(0.5,1.5),
          (5,b,6): P.closed(2,4),
          (5,b,3): P.closed(1,3)
          }
    
    G = fsa(Xt,Et,Tt,X0t,Xmt,Sigobs=sigobst,name="$G_T$")
    GT = tia(G,mut)
    # G_T in bearmer format
    ti_draw(GT)

    # G_T in bearmer format
    ti_draw(GT, 'beamer')

    # G_T in figure format
    ti_draw(GT, 'figure')

    # G_T in figurecolor format
    ti_draw(GT, 'figurecolor')
    """
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
    return


def DP(Gt,state):
    """
    Return all the possible detectable paths that start from the given 
    state of a time-interval automaton.
    
    -------
    Example

    syms('a b u')

    # automaton definition G_T
    Xt = [0, 1, 2, 3, 4, 5, 6]
    Et =[a,b,u]
    sigobst = [a,b]
    X0t = [0]
    Xmt = []
    Tt = [(0,u,1),(0,a,4),(1,a,2),(2,b,3),(4,u,5),(5,b,6),(5,b,3)]
    
    mut = {(0,u,1): P.closed(0,2),
           (0,a,4): P.closed(1,3.5),
           (1,a,2): P.closed(1,3),
           (2,b,3): P.closed(2,3),
           (4,u,5): P.closed(0.5,1.5),
           (5,b,6): P.closed(2,4),
           (5,b,3): P.closed(1,3)
          }
    
    G = fsa(Xt,Et,Tt,X0t,Xmt,Sigobs=sigobst,name="$G_T$")
    GT = tia(G,mut)

    # Print the detectable path from state 4
    print(DP(GT, 4))
    """
    G = Gt[0]
    DP = [] 
    Eobs = list(G.Sigobs)
    E_active = list(G.Gamma(state))

    for event in E_active:
        evact = G.delta(state,event)
        if isinstance(evact,frozenset):
            conj = list(evact)
            for el in conj:
                DP += [[(state,event,el)]]
        else:
            conj = evact
            DP += [[(state,event,conj)]]
        
    aux = DP[:]
    while len(aux) != 0:
        for p in DP:
            last_event = p[len(p) - 1][1]
            last_state = p[len(p) - 1][2]
            if p in aux:
                aux.remove(p)
            if last_event not in Eobs:
                DP.remove(p)
                E_actives_p = G.Gamma(last_state)
                for ev in E_actives_p:
                    if type(G.delta(last_state,ev)) == frozenset:
                        for dest_state in G.delta(last_state,ev):
                            p_actual = p[:]
                            p_new = (last_state, ev, dest_state)
                            p_actual.append(p_new)
                            aux.append(p_actual)
                            DP.append(p_actual)
                    else:
                        p_actual = p[:]
                        p_new = (last_state, ev, G.delta(last_state,ev))
                        p_actual.append(p_new)
                        aux.append(p_actual)
                        DP.append(p_actual)
    return DP


def ti_proj(Gt):
    """
    Compute the projection TIA of a time-interval automaton.

    -------
    Example
    
    syms('u')

    # automaton definition G_T
    Xt = [0, 1, 2, 3]
    Et =[a,b,u]
    sigobst = [a,b]
    X0t = [0]
    Xmt = [ ]
    Tt = [(0,u,1),(1,a,2),(2,b,3)]
    
    mut = {(0,u,1): P.closed(0,2),
            (1,a,2): P.closed(1,3),
            (2,b,3): P.closed(2,3)
          }
    
    G = fsa(Xt,Et,Tt,X0t,Xmt,Sigobs=sigobst,name="$G_T$")
    GT = tia(G,mut)
    ti_draw(GT, 'figure')

    # Generate the projection of G_T
    d = ti_proj(GT)
    ti_draw(d, 'figure')
    """

    ## Computes
    ## Calculate the set of destination states for events
    ## observable or non-observable
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
    
    ## Internal function to generate the observable projection
    ## of the path and the generated Interval
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

    ##  MAIN TI PROJ ##
    G = Gt[0]
    mu_old = Gt[1]
    X0_new = list(G.X0)
    X_new = X0_new[:]
    E_new = list(G.Sigobs)
    T_new = []
    mu_new = dict()
    
    # Adds to X_new the destination states of transitions labeled by observable events
    X_new += SubsetX(G, 'Obs')
    
    # Calculates the Detectable Path for each state in X with Gamma > 0
    path = []
    for x in X_new:
        if list(G.Gamma(x)):
            dp = DP(Gt,x)
            path += dp
            
    # Projects the Detectable Paths, creating the new transitions and their associated I
    proj = []
    for p in path:
        p_proj, proj_mu = ObservableProj(p,mu_old)
        proj.append(p_proj)
        mu_new.update(proj_mu)
    T_new = proj

    #marking states
    Xm_new = []
    for x in X_new:
        if x in G.Xm:
            Xm_new +=[x]
    
    # Creates the projection automaton and combines it with the time function mu_new
    Pi_Gt = fsa(X_new, E_new, T_new, X0_new, Xm_new, name = '$\Pi(' + G.name + ')$')
    Pi_Gt.setgraphic(style='observer')
    Gproj = (Pi_Gt,mu_new)
                    
    return Gproj


## Returns the intersection of all intervals
def mu_intersection(mu, intervals):
    inter = P.closedopen(0,P.inf)
##    print("o que é isso: " + str(intervals))
    for i in intervals:
##        print("isso que esta tentando fazer a operação: " + str(i))
        inter = inter & i
        
    return inter

def ti_equi_det(Gt):
    """
    Remove nondeterminism from a TIA which are:
        There is more than one initial state;

        There is more than one transition originating from a state,
        labeled with the same event, and whose time intervals 
        are not disjoint

    -------
    Example

    Xt = [0, 1, 2, 3, 4, 5, 6]
    Et =[a,b,c]
    sigobst = [a,b,c]
    X0t = [0,1]
    Xmt = [ ]
    Tt = [(0,c,1),(0,c,4),(1,a,2),(2,b,3),(4,c,5),(5,b,6),(5,b,3)]
    
    mut = {(0,c,1): P.closed(0,2),
            (0,c,4): P.closed(1,3.5),
            (1,a,2): P.closed(1,3),
            (2,b,3): P.closed(2,3),
            (4,c,5): P.closed(0.5,1.5),
            (5,b,6): P.closed(2,4),
            (5,b,3): P.closed(1,3)
          }
    
    G = fsa(Xt,Et,Tt,X0t,Xmt,Sigobs=sigobst,name="$G_T$")
    GT = tia(G,mut)
    ti_draw(GT, 'figure')

    #Generate the deterministic equivalent of G_T
    d = ti_equi_det(GT)
    ti_draw(d, 'figure')
    """

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


## Returns a dictionary with disjoint transitions, 
# including new states and their associated I.
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


## Renames observable events at the end of a Detectable Path and
# the states of Glt with the events that lead to it.
## ret = 'GRF': Returns Glt with renamed states (Renamed Final).
## ret = 'GRI': Returns Glt with renamed transitions (Renamed Intermediate).
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

## Under Construction: Unpacks the renamed states of Gdtr, renaming 
# states and creating or replacing transitions
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


def ext_ti_product(self,other):
    G1, Mu1 = self
    G2, Mu2 = other    
    sigp_1 = G1.Sigma - G2.Sigma
    sigp_2 = G2.Sigma - G1.Sigma
    
    for sig in sigp_2:
        for x in G1.X:
            G1 = G1.addselfloop(x,sig)
            Mu1[(x,sig,x)] = P.closedopen(0,P.inf)
    for sig in sigp_1:
        for x in G2.X:
            G2 = G2.addselfloop(x,sig)
            Mu2[(x,sig,x)] = P.closedopen(0,P.inf)

    extprod = ti_product(tia(G1,Mu1),tia(G2,Mu2))
    extprod[0].setgraphic(style='observer')
    return extprod



def ti_product(self, other):
    """
    Performs the product composition between two TIAs, where 
    transitions synchronize on common events, and the intersection
    of their respective time intervals is computed.

    -------
    Example

    from deslab import *

    # automaton definition G1
    Xt1 = [0, 1, 2]
    Et1 =[a,b]
    sigobst1 = [a,b]
    X0t1 = [0]
    Xmt1 = [ ]
    Tt1 = [(0,a,1),(1,b,2)]
    
    mut1 = \{(0,a,1): P.closed(1,3),
            (1,b,2): P.closed(0,2)
        \}
    
    G1 = fsa(Xt1,Et1,Tt1,X0t1,Xmt1,Sigobs=sigobst1,name="$G_{1}$")
    G_tia1 = tia(G1,mut1)
    ti_draw(G_tia1, 'figure')
    
    # automaton definition G2
    Xt2 = [0, 1, 2, 3, 4]
    Et2 =[a,b,c]
    sigobst2 = [a,b,c]
    X0t2 = [0]
    Xmt2 = [ ]
    Tt2 = [(0,a,1),(1,c,2),(2,b,3),(0,c,4)]
    
    mut2 = {(0,a,1): P.closed(2,4),
            (1,c,2): P.closed(1,3),
            (2,b,3): P.closed(1,4),
            (0,c,4): P.closed(0,4)
        }
    
    G2 = fsa(Xt2,Et2,Tt2,X0t2,Xmt2,Sigobs=sigobst2,name="$G_{2}$")
    G_tia2 = tia(G2,mut2)
    ti_draw(G_tia2,'figure')

    # Product
    ti_draw(ti_product(G_tia1,G_tia2),'figure')
    """
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
        p1, p2 = p     
        for sigma in Gamma_prod(p):
            Q1_n, Q2_n = delta1(p1, sigma), delta2(p2,sigma)
            if type(Q1_n) == tuple:
                Q1_n = set([Q1_n])
            if type(Q2_n) == tuple:
                Q2_n = set([Q2_n])
            try:
                type1 = (type(Q1_n) == str or type(Q1_n) == int)
                type2 = (type(Q2_n) == str or type(Q2_n) == int)
                if type1 and not type2:
                    Q = list()
                    for q2 in Q2_n:
                        Q.append((Q1_n,q2))
                elif not type1 and type2:
                    Q = list()
                    for q1 in Q1_n:
                        Q.append((q1,Q2_n))
                elif type1 and type2:
                    Q = [(Q1_n,Q2_n)]
                elif not type1 and not type2:
                    Q = list((q1,q2) for q1 in Q1_n  for q2 in Q2_n)
                
                
                # Q contains the tuples for reached states
                for q in Q:
                    if Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])] != P.empty():
                        if q not in X_p:
                            # updating states of composite automaton and stack
                            X_p = X_p | set([q])           
                            S.append(q)
                            q1, q2 = q
                            # updating marked states
                            if (q1 in Xm1) & (q2 in Xm2):
                                Xm_p = Xm_p | set([q]) 
                        transition.append([p,sigma,q])
                        if not mu_p.get((p,sigma,q)):
                            mu_int = Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])]
                            mu_p[(p,sigma,q)] = mu_int
            except:
                q = (Q1_n,Q2_n)
                if Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])] != P.empty():
                    if q not in X_p:
                        # updating state and stack                                                     
                        X_p = X_p | set([q])           
                        S.append(q)
                        q1, q2 = q
                        if (q1 in Xm1) & (q2 in Xm2):
                            # updating marked states
                            Xm_p =Xm_p | set([q])
                    transition.append([p,sigma,q])
                    if not mu_p.get((p,sigma,q)):
                        mu_int = Mu1[(p1,sigma,q[0])] & Mu2[(p2,sigma,q[1])]
                        mu_p[(p,sigma,q)] = mu_int

    # finally, we build the composite automaton
    prodnondet = fsa(list(X_p), Sigma_p, transition, X0_p, list(Xm_p), table=table_p,
                  Sigobs = Sigobs_p, Sigcon=Sigcon_p,name='untitled')
    tiprod = tia(prodnondet,mu_p)
    return tiprod

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

## Same function of diagnosis toolbox, but for TIA
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


def ti_complement(Gt):
    """
    Generates the TIA that marks the complement of the marked 
    time-interval language of a given TIA.

    -------
    Example
    
    # automaton definition G_T
    Xt1 = [0, 1]
    Et1 =[a]
    sigobst1 = [a]
    X0t1 = [0]
    Xmt1 = [ ]
    Tt1 = [(0,a,1)]
    
    mut1 = {(0,a,1): P.closed(1,3)
           }
    
    G1 = fsa(Xt1,Et1,Tt1,X0t1,Xmt1,Sigobs=sigobst1,name="$G_{1}$")
    G_tia1 = tia(G1,mut1)
    ti_draw(G_tia1,'figure')

    # Complement
    ti_draw(ti_complement(G_tia1),'figure')
    """
    # Unmark all states and add the dump state xd
    G = Gt[0]
    mu_new = Gt[1].copy()
    NewMarked = list(G.X - G.Xm) + ['dump']
    G = G.addstate('dump')
    G = G.setpar(Xm=NewMarked)
    
    Gtcomp = tia(G,mu_new)

    for x in G.X:
        trans = []
        active = []
        # Set of transitions that leave x
        for t in G.transitions():
            if t[0]==x:
                trans += [t]
                active += [t[1]] # Active events
        # Compute the complement of the union of the intervals associated 
        # with the event ev that leaves the state x
        for ev in active:
            union = P.empty()
            for t in trans:
                if t[1] == ev:
                    union = union | mu_new[t]
            comp = P.closedopen(0,P.inf) - union # Complementary interval
            # Create a transition with the event ev and the complementary interval
            # comp that leads to the dump state
            newt = (x,ev,'dump') #new transition
            mu_new[newt] = comp # Interval associated with the new transition
            G = G.addtransition(newt)
            Gtcomp = tia(G, mu_new)

        #Transitions to the dump state with events not defined in the state
        for ev in list(set(G.Sigma) - set(active)):
            newt = (x,ev,'dump')
            mu_new[newt] = P.closedopen(0,P.inf)
            G = G.addtransition(newt)
            Gtcomp = tia(G,mu_new)

    # Self-loops in the dump state
    for ev in list(G.Sigma):
        newt = ('dump',ev,'dump')
        mu_new[newt] = P.closedopen(0,P.inf)
        G = G.addtransition(newt)
        Gtcomp = tia(G,mu_new)
        

    return Gtcomp
            




    
