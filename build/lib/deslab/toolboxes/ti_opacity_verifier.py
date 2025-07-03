from deslab import *


## Performs the product operation between two TIA while adding labels
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

## Performs the product operation between two TIA while adding labels Grev
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

def TLBO(Gst,Gnst):
    """
    Verify the property of timed language-based opacity (TLBO) in TIA

    -------
    Example

    syms('0 1 2 3 4 a b')


    # automaton definition G_ST
    Xs = [0,1,2,3,4]
    sigobss = [a,b]
    Es = [a,b]
    X0s = [0]
    Xms = [2,3,4]
    Ts =[(0,a,1),(1,b,2),(0,a,3),(0,b,4)]
    
    mus = {
         (0,a,1) : P.closed(0,2),
         (1,b,2) : P.closed(1,5),
         (0,a,3) : P.closed(0,1),
         (0,b,4) : P.closed(2,4)
         }
    
    Gs = fsa(Xs,Es,Ts,X0s,Xms,Sigobs = sigobss,name="$G_s$")
    Gst = tia(Gs,mus)
    
    # automaton definition G_NST
    Xns = [0,1,2,3,4]
    sigobsns = [a,b]
    Ens = [a,b]
    X0ns = [0]
    Xmns = [2,3,4]
    Tns =[(0,a,1),(1,b,2),(0,a,3),(0,b,4)]
    
    muns = {
         (0,a,1) : P.closed(0,5),
         (1,b,2) : P.closed(1,4),
         (0,a,3) : P.closed(0,2),
         (0,b,4) : P.closed(0,1)
         }
    
    Gns = fsa(Xns,Ens,Tns,X0ns,Xmns,Sigobs = sigobsns,name='$G_{ns}$')
    Gnst = tia(Gns,muns)
    
    # Print the language based opacity property
    print(TLBO(Gst,Gnst))    
    """
    #projection of Gst
    Gstproj = ti_proj(Gst)
    Gstproj[0].name = 'Proj(Gst)'


    #projection of Gnst 
    Gnstproj = ti_proj(Gnst)
    Gnstproj[0].name = 'Proj(Gnst)'

    #observer of Gnst 
    Gnsto = ti_equi_det(Gnstproj)
    Gnsto[0].name = 'Obs(Gnst)'
    

    #complement of Obs(Gnst) 
    Gnstocomp = ti_complement(Gnsto)
    Gnstocomp[0].name = 'Comp(Obs(Gnst))'

    #G Obfuscated 
    Gobf = ti_label_obf(Gstproj,Gnstproj)

    Gobf = tia(coac(Gobf[0]),Gobf[1])

    #G Revealed
    Grev = ti_label_rev(Gstproj,Gnsto)

    Grev = tia(coac(Grev[0]),Grev[1])


    if isitempty(Gobf[0]):
        decision = 'Not opaque'
        return decision
        sys.exit()
    
    if isitempty(Grev[0]):
        decision = 'TISO'
        return decision
        sys.exit()


    #PART 1 VERIFICATION
    #Generate the verifier automaton with labels on the transitions
    Gv,dict_v, dict_state = verifierTLBO(Gst[0],Gnst[0],Gobf,Grev)

    Gv_label = tia(Gv,dict_v)
    dict_original = deepcopy(Gv_label[1])
    Gv_label_original = tia(Gv_label[0],dict_original)
    
    #PART 2 VERIFICATION
    #Assign a label to each transition of the verifier
    #Following the rules that define this

    for t in Gv.transitions():
        labels = dict_v[t]
        l_obf = labels[0]
        l_rev = labels[1]

        #1
        if 'empty' in l_rev: 
            dict_v[t] = 'co'
            
        #2    
        if 'empty' in l_obf: 
            dict_v[t] = 'cr'  

        
        if ('empty' not in l_obf) and (('pr' in l_rev) or ('cr' in l_rev)): #and ('co' not in l_rev) and ('po' not in l_rev):
            dict_v[t] = 'por'

        if ('empty' not in l_obf) and (('po' in l_rev) or ('co' in l_rev)) and ('cr' not in l_rev) and ('pr' not in l_rev):
            dict_v[t] = 'co'



    #PART 3
    # Verification decision
    flag_obf = 0
    flag_part = 0
    flag_rev = 0

    Gobf_strings = Gv
    Gpart_strings = Gv

    for t in Gv.transitions():
        if dict_v[t] == 'cr':
            flag_rev = 1
        
        if (dict_v[t] != 'co'):
            Gobf_strings = Gobf_strings.deletetransition(t)

        if (dict_v[t] == 'cr'):
            Gpart_strings = Gpart_strings.deletetransition(t)

    if not isitempty(trim(Gobf_strings)):
        flag_obf = 1

    for t in Gpart_strings.transitions():
        if dict_v[t] == 'por':
            flag_part = 1

    if (flag_part == 1)and (flag_rev == 0):
        decision =  'TDSO'

    if (flag_obf == 1) and (flag_rev == 1):
        decision = 'TIWO'

    if (flag_obf == 0) and (flag_part == 1) and (flag_rev == 1):
        decision = 'TDWO'

    return decision

