from deslab import *

def current_state_op(G, Xs, Xns=[]):
    """
    Return True if the current state opacity exists, otherwise return False.
    Xs represents the secret states, and Xns represents the non-secret states.
    
    -------
    Example
    syms('q0 q1 q2 q3 q4 a1 b1 c1 d1')
    table = [(a1,'a'),(b1,'b'),(c1,'c'),(d1,'d'),
             (q1,'q_1'),(q2,'q_2'),(q3,'q_3'),(q0,'q_0'),(q4,'q_4')]
    X = [q0,q1,q2,q3,q4]
    Sigma = [a1,b1,c1,d1]
    Sigmao = [a1,c1,d1]
    X0 = [q0]
    Xm = []
    T =[(q0,a1,q1),(q1,b1,q2),(q1,d1,q3),(q2,c1,q2),(q2,d1,q4), (q3,b1,q4)]
    G = fsa(X,Sigma,T,X0,Xm,table,Sigmao, name='$G$') 
    
    xs = [q3]
    xns = [q4]
    is_current_state_opaque = current_state_op(G, xs, xns)
    """
    
    SigmaOb = list(G.Sigobs)
    
    G2 = observer(G, SigmaOb)
    
    
    all_states = G2.X

    if Xns==[]:
        Xns = [x for x in G.X if x not in Xs]
    
    for  states in all_states:
        secret_state = False
        nonsecret_state = False
        for x in states:
            if x in Xs:
                secret_state = True
            if x in Xns:
                nonsecret_state = True
        if secret_state and not nonsecret_state:
            return False
    return True

def inverse_automaton(G):
    trans = transitions(G)
    for x1, t1, x2 in trans:
        G = G.deletetransition((x1, t1, x2))
        G = G.addtransition([x2,t1,x1])

    return G

def initial_state_opac(G, Xs, Xns=[]):
    """
    Return True if the initial state opacity exists, otherwise return False.
    Xs represents the secret states, and Xns represents the non-secret states.
    
    -------
    Example
    syms('q0 q1 q2 q3 q4 a1 b1 c1 d1')
    table = [(a1,'a'),(b1,'b'),(c1,'c'),(d1,'d'),
             (q1,'q_1'),(q2,'q_2'),(q3,'q_3'),(q0,'q_0'),(q4,'q_4')]
    X = [q0,q1,q2,q3,q4]
    Sigma = [a1,b1,c1,d1]
    Sigmao = [a1,b1,c1]
    X0 = [q0,q1,q2]
    Xm = []
    T =[(q0,a1,q1),(q1,b1,q2),(q1,d1,q3),(q2,c1,q2),(q2,d1,q4), (q3,b1,q4)]
    G = fsa(X,Sigma,T,X0,Xm,table,Sigmao, name='$G$') 
    
    xs = [q2]
    xns = [q0,q1,q3,q4]
    is_initial_state_opaque = initial_state_opac(G, xs, xns)
    """
    
    SigmaOb = list(G.Sigobs)
    
    if Xns==[]:
        Xns = [x for x in G.X if x not in Xs]


    if sorted(G.X) == sorted(G.X0):
        allow_states = list(G.X0)
        G1 = G.setpar(Xm = G.X)
        G2 = inverse_automaton(G)
        G3 = observer(G2, SigmaOb)
    else:
        if not any(x in list(G.X0) for x in Xs):
            return False
        else:
            allow_states = list(G.X0)
            G1 = G.setpar(X0 = G.X, Xm = G.X0)
            G2 = inverse_automaton(G1)
            G3 = observer(G2, SigmaOb)


    marked_states = G3.Xm
    for states in marked_states:
        secret_state = False
        nonsecret_state = False
        for x in states:
            if (x in Xs) and (x in allow_states):
                secret_state = True
            if (x in Xns) and (x in allow_states):
                nonsecret_state = True
        if secret_state and not nonsecret_state:
            return False
    return True
    
def language_based_opac(G1, G2, SigmaO):
    """
    Return True if the language based opacity exists, otherwise return False.
    G1 represents the secret language, G2 represents the non-secret language and 
    SigmaO represents the observable events.
    
    ------------
    Example
    syms('q0 q1 q2 q3 q4 a1 b1 c1 d1')
    table = [(a1,'a_1'),(b1,'b_1'),(c1,'c_1'),(d1,'d_1'),(q1,'q_1'),(q2,'q_2'),(q3,'q_3'),(q0,'q_0'),(q4,'q_4')]
    X = [q0,q1,q2,q4]
    Sigma = [a1,b1,d1]
    SigmaO = [a1,d1]
    X0 = [q0]
    Xm = [q4]
    T =[(q0,a1,q1),(q1,b1,q2),(q2,d1,q4)]
    G = fsa(X,Sigma,T,X0,Xm,table, SigmaO,name='$G$')
    
    table1 = [(a1,'a_1'),(b1,'b_1'),(c1,'c_1'),(d1,'d_1'),(q1,'q_1'),(q2,'q_2'),(q3,'q_3'),(q0,'q_0'),(q4,'q_4'),(q5,'q_5')]
    X1 = [q0,q1,q2,q3,q4]
    Sigma1 = [a1,b1,c1,d1]
    X01 = [q0]
    Xm1 = [q4]
    T1 =[(q0,a1,q1),(q1,b1,q2),(q1,d1,q3),(q2,c1,q2),(q2,d1,q4), (q3,b1,q4)]
    G1 = fsa(X1,Sigma1,T1,X01,Xm1,table1,name='$G_1$')
    
    sigma_o = [a1, c1, d1]
    is_language_based_opaque = language_based_opac(G, G1, sigma_o)
    """
    
    SigmaOb = SigmaO
    
    G1o = observer(G1, SigmaOb)
    G2o = observer(G2, SigmaOb)
    G3 = coac(product(G1o, G2o))
    G1oc = complement(G1o)
    G4 = coac(product(G3, G1oc))

    #draw(G1, G2)
    #draw(G1o,G2o)
    #draw(G3, G1oc, G4)
    
    if G3.X != []:
        if langdiff(G4,G1o).X != []:
            return True
    else:
        return False
    
def initial_final_state_opac(G, xsp, xnsp=[]):
    """
    Return True if the initial-final state opacity exists, otherwise return False.
    xsp represents the secret states pairs, and xnsp represents the non-secret states pairs.
    
    ------------
    Example
    syms('q0 q1 q2 q3 q4 a1 b1 c1 d1 e1')
    table = [(a1,'a_1'),(b1,'b_1'),(c1,'c_1'),(d1,'d_1'),(q1,'q_1'),(q2,'q_2'),(q3,'q_3'),(q0,'q_0'),(q4,'q_4')]
    X = [q0,q1,q2,q3]
    Sigma = [a1,b1,e1]
    SigmaO = [a1,b1]
    X0 = [q0,q2]
    Xm = []
    T =[(q0,a1,q0),(q0,e1,q2),(q1,b1,q0),(q2,a1,q1),(q1,e1,q3),(q3,b1,q1)]
    G = fsa(X,Sigma,T,X0,Xm,table, SigmaO,name='$G$')
    
    
    xps = [('q2','q1')]
    xnps = [('q2','q2')] 
    print(initial_final_state_opac(G, xps, xnps))
    """
    
    m0 = []
    keys = []
    m_aux = []
    m_aux1 = []
    SigmaNob = []
    trans_nob = []
    remove_keys = []
    X0 = G.X0
    trans = list(G.transitions())
    contador = 0
    states = {}
    states_i = {}
    states_aux = {}

    SigmaOb = list(G.Sigobs)


    if xnsp==[]:
        xnsp_aux = [(a, b) for a in G.X for b in G.X]
        xnsp = [xp for xp in xnsp_aux if xp not in xsp]

    # building m0
    for i in X0:
        aux = (i,i)
        m0.append(aux)
    

    for begin, event, end in trans:
        if event not in SigmaOb:
            SigmaNob.append(event)
        if begin in X0 and event not in SigmaOb:
            aux = (begin, end)
            m0.append(aux)

    nome_lista = f'm_{contador}'
    states_i[nome_lista] = m0

    for begin, event, end in trans:
        if event in SigmaNob:
            aux = (begin, end)
            trans_nob.append(aux)

    count = 0
    # building others m
    while len(states_i) != 0:
        count = count + 1
        keys = []
        remove_keys = []

        for i in states_i.keys():
            remove_keys.append(i)
            states[i] = states_i[i]
            
            for e in SigmaOb:
                m_aux = []
                m_aux1 = []
                m_existe = False
               
                for state in states_i[i]:
                    for begin, event, end in trans:
                        if begin == state[1] and e == event:
                            aux = (state[0],end)
                            m_aux.append(aux)

                while len(m_aux) != 0:
                    m_aux2 = m_aux
                    m_aux = []
                    for dupla in m_aux2:
                        m_aux1.append(dupla)
                        for begin,end in trans_nob:
                            if dupla[1] == begin:
                                aux = (dupla[0],end)
                                m_aux.append(aux)
             
                for nome_lista, lista in states.items():
                    if m_aux1 == lista:
                        m_existe = True
                
                if m_existe == False:
                    contador = contador+1
                    nome_lista = f'm_{contador}'
                    states_aux[nome_lista] = m_aux1

        for i in states_aux.keys():
            keys.append(i)
        
        for i in keys:
            states_i[i] = states_aux.pop(i)

        for key in remove_keys:
            del states_i[key]
     
    
    for i in states.keys():
        estado_s = False
        estado_ns = False

        for par in states[i]:
            if par in xsp:
                estado_s = True
            if par in xnsp:
                estado_ns = True

        if estado_s and not estado_ns:
            return False
    if not (estado_s and not estado_ns):
        return True
        
