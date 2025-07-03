from deslab import *
from networkx import selfloop_edges

syms('q0 q1 q2 q3 q4 q5 a1 b1 c1 d1 e1 a b c d e f x y t1 t2 su')
#s = list(s) retirar frozenset
#deepcopy para salvar o estado original e poder alterar o automato
#lexgraph alphamap(G1) pra chregar
#olhar variaveis
#os estados em Y ou Z tem que ser ajustados, olhar o final do que tem no artigo




def verifier_estimator(G, Xs):
    events = G.Sigma
    G = G.addevent('new_event')
    G = G.setpar(Sigobs=events)
    e = observer(G, G.Sigobs)

    #montando E_d
    e_d = e
    all_states = e.X
    for  states in all_states:
        secret_state = False
        nonsecret_state = False
        for x in states:
            if x in Xs:
                secret_state = True
            else:
                nonsecret_state = True
        if secret_state and not nonsecret_state:
            e_d = e_d.deletestate(states)
    e_d = ac(e_d)

    #montando E_f
    trans = e.transitions()
    e_f = e

    for x1,ev,x2 in trans:  
        if (x1 != x2) and (ev in G.Sigobs):
            e_f = e_f.addtransition([x1,f,x2])
    for i in e_f.Sigma:
        evi = i + "i"
        e_f = e_f.addevent(evi)
    for states in all_states:
        s_events = []
        for x1,ev,x2 in trans:  
            if states == x1 and x1 == x2:
                s_events.append(ev)
        for ev in e.Sigma:
            if ev not in s_events:
                evi = ev + "i"
                e_f = e_f.addtransition([states,evi,states])

    return (e,e_d,e_f)

def verifier_parallel_composition(e,e_d, e_f):
    
    new_states = []
    all_states = []
    initial_state = []
    v_fvo = []
    v_fvi = []
    v_fve = []

    x01 = next(iter(e_d.X0))
    x02 = next(iter(e_f.X0))
    
    trans1 = e_d.transitions()
    trans2 = e_f.transitions()

    aux = (list(x01)[0],list(x02)[0])
    new_states.append(aux)
    initial_state.append(aux)
    new_xv = []
    while len(new_states) != 0:
        
        

        #isso pode dar errado por causa do frozenset
        for state in new_states:
            state_d = state[0]
            state_f = state[1]
            all_states.append(state)
            
            for i in e_d.Sigma:
                #transicao em e_d
                trans_d = [group for group in trans1 if state_d in group[0] and group[1] == i]

                #transicao normal
                trans_fo = [group for group in trans2 if state_f in group[0] and group[1] == i] 
                #transicao self loop criado
                trans_fi = [group for group in trans2 if state_f in group[0] and group[1] == i+"i"]
                
                #transicao epsilon - pega a mesma info varias vezes mas ta mais bonito deixar aqui
                trans_fe = [group for group in trans2 if state_f in group[0] and group[1] == "f"]
                
                
                if trans_d != [] and trans_fo != []:
                    x1 = (state_d,state_f)
                    x2 = (next(iter(trans_d[0][2])),next(iter(trans_fo[0][2])))
                    v_fvo.append([x1,i,x2])
                    new_xv.append(x2)
                if trans_d != [] and trans_fi != []:
                    x1 = (state_d,state_f)
                    x2 = (next(iter(trans_d[0][2])),state_f)
                    v_fvi.append([x1,i,x2])
                    new_xv.append(x2)
                if trans_fe != []:
                    x1 = (state_d,state_f)
                    x2 = (state_d,next(iter(trans_fe[0][2])))
                    v_fve.append([x1,i,x2])
                    new_xv.append(x2)
        
        aux1 = set()
        new_xv_aux = []
        for sublist in new_xv:
            if sublist not in new_xv_aux:
                new_xv_aux.append(sublist)

        new = []
        new_xv_aux = new_xv_aux
        for item in new_xv_aux:
            if item not in all_states:
                new.append(item)

        new_states = []
        new_states = new


    T = []


    for sublist in v_fvo:
        T.append(sublist)
    for sublist in v_fvi:
        T.append(sublist)
    for sublist in v_fve:
        T.append(sublist)

    v = fsa(all_states,list(e.Sigma),T,initial_state,[],name='$V$')
    #draw(v)
    return (v,v_fvo,v_fvi,v_fve)

def unfolded_verifier(v,fvo,fvi,fve):
    if True:
        y = []
        z = []
        final_y = []
        final_z = []
        f_zy = []
        f_yz = []
        vu_sigma = []
        states_v = v.X
        sigma_v = v.Sigma
        
        y0 = v.X0
        y = y0
    
    while y != []:
        for state_y in y:
            
            for i in sigma_v:    
                
                v_fvo_e = [group for group in fvo if group[0] == state_y and group[1] == i]
                v_fve_e = [group for group in fve if group[0] == state_y and group[1] == i]
                
                if v_fvo_e != []:
                    f_yz.append([state_y,i,(state_y,i)])
                    if (state_y,i)not in z:
                        z = z + [(state_y,i)] 
                elif v_fve_e != []:
                    f_yz.append([state_y,i,(state_y,i)])
                    if (state_y,i)not in z:
                        z = z + [(state_y,i)]

        final_y += y
        y = []
        v_aux = fsa()
        auto_vazio = True
        trans = {}
        #state_z[0] eh o estado e state_z[1] o evento
        if True:
            for state_z in z:
                xv_aux = []
                vz = v
                vz = vz.setpar(X0 = state_z[0])
                vz_ac = ac(vz)
                t_vz = vz_ac.transitions() #acho que pode dar problema                
                aux2 = []
                for x1,ev,x2 in t_vz:
                    #print("x1",x1)
                    if ev == state_z[1]:      
                        for group in fvo:
                            if group[0] == x1 and group[1] == ev and group[2] == x2:
                                xv_aux.append([group[0],group[2]])

                #problema eh aqui
                for state,final in xv_aux:
                    if True: #state_z[0] != state:
                        v_aux = fsa(list(states_v),list(sigma_v),fvi,state_z[0],[state])
                        v_aux = coac(v_aux)
                        auto_vazio = v_aux.empty
                        xv_linha1 = state
                        #final_linha = final
                        if not auto_vazio:
                            trans = lexgraph_alphamap(v_aux)
                            aux2.append([state,final,trans])
                    
                if xv_aux != [] and aux2 != []:
                    for xv_linha,final_linha,trans in aux2:
                        #trans = lexgraph_alphamap(v_aux)
                        xv_way = trans[xv_linha]
                        vu_sigma.append(xv_way)
                        f_zy.append([state_z,xv_way,final_linha])
                        if final_linha not in final_y:
                            y.append(final_linha)
                

                for group in fve:
                    if group[0] == state_z[0] and group[2] == xv_linha1:
                        if xv_linha not in final_y:
                            vu_sigma.append('f')
                            f_zy.append([state_z,'f',xv_linha1])
                            y.append(xv_linha)
                
        #print("olhar aqui Y\n", y)
        #print("f_zy\n",f_zy)

        final_z += z
        z = []

    #print(vu_sigma)
    #print("Y\n", final_y + final_z)
    #print("Z\n", final_z)
    #print("F_zy \n", f_zy+f_yz)
    #print("F_yz\n",f_yz)

    return final_y,final_z,f_zy,f_yz, list(set(vu_sigma))
    
def all_edit_structure_c(vu, edit_const):
    
    trans_old = []
    vum = vu
    for state in edit_const:
        vum = vum.deletestate(state)
    
    vum = ac(vum)
    
    all_states = vum.X
    trans_new = vum.transitions()
    state_active = []
    
    for x1,ev,x2 in trans_new:
        if x1 not in state_active:
            state_active.append(x1)
    
    while trans_old != trans_new:
        trans_old = trans_new
        
        for state in all_states:
            if state not in state_active:
                vum = vum.deletestate(state)

        vum = ac(vum)
        
        all_states = vum.X
        trans_new = vum.transitions()
        state_active = []
        for x1,ev,x2 in trans_new:
            if x1 not in state_active:
                state_active.append(x1)
        
    return vum

def string_run_and_projection(aes_c,y):
    all_states = aes_c.X
    y0 = list(aes_c.X0)
    
    aes_t = aes_c
    
    for state in all_states:
        state_real = False
        for s in y:
            if state == s:
                state_real = True
        if not state_real:
            aes_t = aes_t.deletestate(state)
        
    aes_t = ac(aes_t)
    
    trans = aes_t.transitions()
    
    current_state = y0[0]   
    s_run = []
    s_run_p = []
    while trans != []:
        
        for t in trans:
            if t[0] == current_state:
                initial_state = t[0]
                e_aux = t[1]
                current_state = t[2]
                x = t
                break
        
        trans.remove(x)
        for t in trans:
            if t[0] == current_state:
                gamma_aux = t[1]
                current_state = t[2]
                x = t
                break
        
        trans.remove(x)
        
        
        if gamma_aux != 'epsilon':
            s_run.append([initial_state,(gamma_aux+e_aux),current_state])
            s_run_p.append([initial_state,(e_aux),current_state])
        else:
            s_run.append([initial_state,(e_aux),current_state])
            s_run_p.append([initial_state,(e_aux),current_state])
            
    #print("s_run ", s_run)
    #print("s_run_p ", s_run_p)
    
    final_transition = []
    for t in s_run:
        for p in s_run_p:
            if t[0]==p[0] and t[2]==p[2]:
                final_transition.append((t[0],(t[1],p[1]),(t[2],(t[1],p[1]))))
    
    ordem = 0
    current_state = y0[0]
    trans_ordenadas = []
    while final_transition != []:
        if trans_ordenadas == []:
            for t in final_transition:
                if t[0] == current_state:
                    trans_ordenadas.append(t)
                    x = t
                    current_state = t[2]
                    break
        else:
            for t in final_transition:
                if t[0] == current_state[0]:
                    new_t2 = (t[2][0],(current_state[1][0]+t[2][1][0],current_state[1][1]+t[2][1][1]))
                    trans_ordenadas.append([t[0],t[1],new_t2])
                    
                    x = t
                    current_state = new_t2
                    break
        final_transition.remove(x)
    
    return trans_ordenadas
    
def all_edit_structure_t(aes_c):
    if True:
        y0 = list(aes_c.X0)
        trans = aes_c.transitions()
        way_dic = {y0[0] : [y0[0]]}#rever isso se tiver mais de um estado inicial
        aux_dic = {}
        way_end = {}
    
    while way_dic != {}:
        
        if True:
            del_key = []
            aux_dic = {}
            trans_def = False
        for key, state in way_dic.items():
            for x in trans:
                if key == x[0]:
                    aux_dic[x[2]] = state +  [x[0]]
                    trans_def = True
                    #states_z.append(x[2])
            if not trans_def:
                way_end[len(way_end)+1] = state +  [key]
            del_key.append(key)
        
        for key in del_key:
            way_dic.pop(key)
         
        del_key = []
        for key, state in aux_dic.items():
            for x in trans:
                if key == x[0]:
                    way_dic[x[2]] = state +  [x[0]]
            del_key.append(key)
        
        for key in del_key:
            aux_dic.pop(key)
        
        del_key = []
        for key, state in way_dic.items():
            if key in state:
                way_end[len(way_end)+1] = state +  [key]
                del_key.append(key)
        
        for key in del_key:
            way_dic.pop(key)
    
    for key, state in way_end.items():
        state.pop(0)
        way_end[key] = state
    
    new_y = []
    for key,state in way_end.items():
        new_y.append(string_run_and_projection(aes_c,state))
    
    #print("new_y", new_y)
    fw_dic = {}
    n_way = 0
    if len(new_y)>1:        
        for way in new_y:
            fw_dic[n_way] = []
            hold = False
            for item in way:
                if not hold:
                    initial = item[0]
                final = item[2]
                trans = item[1][1]
                trans2 = item[1][0]
                trans2 = trans2.replace(trans,'')
                if trans2 == '':
                    trans2 = 'epsilon'
                if not hold:
                    aux = [initial,trans,(initial,trans)]
                else:
                    aux = [initial,trans,(initial[0],trans)]
                #print("aux  1 ",aux)
                fw_dic[n_way] = fw_dic[n_way] + [aux]
                
                
                if not hold:
                    aux = [(initial,trans),trans2,final]
                else:
                    aux = [(initial[0],trans),trans2,final]
                #print("aux  2 ",aux)
                fw_dic[n_way] = fw_dic[n_way] + [aux]
                
                initial = item[2]
                hold = True

            n_way += 1
    else:
        fw_dic[0]=[]
        for item in new_y:
            initial = item[0]
            final = item[2]
            trans = item[1][1]
            trans2 = item[1][0]
            trans2 = trans2.replace(trans,'')
            
            aux = [initial,trans,(initial,trans)]
            fw_dic[0] = fw_dic[0] + [aux]
            aux = [(initial,trans),trans2,final]
            fw_dic[0] = fw_dic[0] + [aux]
    
    return fw_dic #fazer o automato tbm
     
def verifier_edit_function(runs,xs, x0):
    secret_run = False
    cj = []
    ck = []
    cj_final = []
    ck_final = []
    leaf = []
    map_runs = []
    used_leaf = []
    aux = []
    aux_run = []
    
    for num in runs:
        run = runs[num]
        secret_run = False
        for state in run:
            if isinstance(state[0][0],tuple):
                if state[0][0][1] in xs:
                    secret_run = True
    
        if secret_run:
            aux = 0
            for state in run:
                if isinstance(state[0][0],tuple):
                    if state[0][0][1] in xs:
                        if aux:
                            ck.append(aux[0][0])
                    else:
                        aux = state
            cj.append(run[-1][-1])
        else:
            ck.append(run[-1][-1])
        
    for j in range(len(cj)):
        for k in range(len(ck)):
            if isinstance(cj[j][1],tuple) and isinstance(ck[k][1],tuple):
                if len(cj[j][1][0])<=len(ck[k][1][0]):
                    ck_aux = ck[k][1][0][-len(cj[j][1][0]):]
                    if ck_aux == cj[j][1][0]:
                        cj_final.append(cj[j])
                        ck_final.append(ck[k])
            
    #remove duplicates cj_final e ck_final
    cj_final = list(set(cj_final))
    ck_final = list(set(ck_final))
    
    leaf = cj_final + ck_final
    
    if leaf:
        for state in leaf:
            for num in runs:
                run = runs[num]
                if run[-1][-1]==state:
                    map_runs.append(run)
                    used_leaf.append(state)
    else:
        print("No edit function")
        return fsa()
                    
    for state in used_leaf:
        leaf.remove(state)
    
    if leaf:
        for state in leaf:
            for num in runs:
                aux = []
                run = runs[num]
                for i in run:
                    if state==i[2]:
                        map_runs.append(run)
                        
                    else:
                        aux.append(run)
                        
    
    state = []
    trans = []
    table = []
    # Percorre a tabela e extrai estados e transições
    for sublista in map_runs:
        table += sublista
        for tupla in sublista:
            state.append(tupla[0])
            trans.append(tupla[1])
            state.append(tupla[2])

    # Remove duplicatas nas listas
    state = list(set(state))
    transicoes = list(set(trans))

    #montar automato
    edit_function_auto = fsa(state,transicoes, table, x0)
    edit_function_auto.setgraphic(style = 'rectangle')
    
    return edit_function_auto

#########################################
def edit_function(G,xs, constrantes = []):
    """
    This function returns an automaton that inserts events to 
    enforce current state opacity to an automaton G, where xs are
    the secret states and constrantes are not allowed states combinations.
    If it returns an empty automaton means the enforcement can not be done.
    
    -------
    Example
    syms('q0 q1 q2 q3 q4 a b c d e')
    table = [(a,'a'),(b,'b'),(c,'c'),(d,'d'),(e,'e'),(q1,'q_1'),(q4,'q_4'),(q5,'q_5'),(q2,'q_2'),(q3,'q_3'),(q0,'q_0')]
    X = [q0,q1,q2,q3,q4,q5]
    Sigma = [a,b,c,d]
    Sigmao = [a,b,c,d]
    X0 = [q0]
    Xm = []
    T =[(q0,d,q1),(q0,a,q4),(q0,b,q5),(q1,a,q2),(q2,b,q3),(q3,c,q0),(q4,b,q5),(q5,c,q0)]
    G = fsa(X,Sigma,T,X0,Xm,table,Sigobs=Sigmao, name='$G$')
    
    x_secret = [q5]
    pp = edit_function(G,x_secret)
     
    """

    
    e,ed,ef = verifier_estimator(G, xs)

    v,fvo,fvi,fve = verifier_parallel_composition(e,ed,ef)

    #a = Y, b = Z
    a,b,c,d,e = unfolded_verifier(v,fvo,fvi,fve)

    vu = fsa(a+b,list(v.Sigma)+e,c+d,list(v.X0),[],name='$Vu$')
    #vu.setgraphic(style = 'rectangle')

    aes_c = all_edit_structure_c(vu,constrantes)
    runs= all_edit_structure_t(aes_c)
    pp_enforce = verifier_edit_function(runs,xs,v.X0)     
    
    return pp_enforce
    
######################################################################################################################

def autoD(auto,OC,Elo=[],loi=0):
    """
    EXAMPLE:

    from DELAYwithSTEPS import *

    syms('s0 s1 s2 s3 s4 s5 s6 a b u c n')
    G_X = [ s1, s2, s3, s4 ]
    G_E = [ u, c, a]
    G_T = [(s1,a,s1),( s1 , u, s2 ),( s2 , c, s3 ),( s1 , c, s4 )]
    G_X0 = [ s1 ]
    G_Xm = [s1, s2, s3, s4]
    G_Econ = [ u, c ] 
    G_Eobs = [ a, c] 
    table = [(s1,'x_1'),(s2,'x_2'),(s3,'x_3'),(s4,'x_4'),(u,r'\mu'),(c,r'\gamma')]
    G = fsa( G_X , G_E , G_T , G_X0 , G_Xm , table = table, Sigcon = G_Econ , Sigobs = G_Eobs ,name ='$G$')

    OC = [(1,[a]) , (0,[c])]
    Elo = [b]

    D = autoD(G,OC,Elo,loi=2)

    draw(G,D)
    
    """

    from collections import deque

    def create_newstate():
        """
        this function can be used to create a new state
        """
        
        countInstance = compCount()
        count = countInstance.counter
        return 'x_{d_{' + str(count) + '}}' # new state

    def rename_event(event,p,c,table,loi=0,ST=False):
        def evname(s):
            if p in s:
                ind = s.index(p) + 1
                if loi != 0:
                    label = s[:ind] + '{' + s[ind:] + p + '{' + '{' + c + '}' + str(loi) + '}' + '}'
                else:
                    label = s[:ind] + '{' + s[ind:] + p + '{' + c + '}' + '}'
            else:
                if loi != 0:
                    label = s + p + '{' + '{' + c + '}' + str(loi) + '}'
                else:
                    label = s + p + '{' + c + '}'
            return label

        s = str(event)
        label = evname(s)
        if ST==False:
            return label
        else:
            s = table.get(event)
            if s == None:
                s = str(event)
            return label,evname(s)

    def create_CHdicts(XC):
        """
        Example:
        XC = [(0.9,[a]) , (0.3,[b,c])]
        dt1,dt2 = create_CHdicts(XC)
        print dt1.items()
        print dt2.items()
        """
        maxDelay = {}
        eventSet = {}
        temp = frozenset([])
        if not XC:
            return maxDelay,eventSet
        for obschannel in XC:
            delay = obschannel[0]
            eset = frozenset(obschannel[1])
            if bool(eset & temp): #if we use List: any(i in eset for i in temp):
                print('ERROR: channel sets are NOT disjoint')
                return {},{}
            temp = temp | eset
            for event in eset:
                maxDelay[event] = delay
                eventSet[event] = eset
        return maxDelay,eventSet
    
    def latexname(q):
        """
        This function takes a state (list) and generates a latex pretty name.
        p = q where q is a queue composed with events  of Sigobs or \nu
        """
        name = ''
        for qi in q:
            if qi in table:
                qi_name = table[qi]
            else:
                qi_name = str(qi)
            name = name + qi_name + ' '
        return name
    
    def newstate(x):
        newState = create_newstate()
        Xlist.update( { newState : x } )
        table.update( { newState : latexname(x) } )
        return newState
     
    def cut(q):
        index = 0
        for qi in q:
            if qi in Sigobs:
                return q[index:]
            index = index+1
        return [_nu_]
    
    def rep(q,i):
        p = q.copy()
        k = len(p)
        if i<0 or i>(k-1):
            print('Error using function rem')
        else:
            p[i] = _nu_
        return p
            

    # Step 0:
    Sigma = auto.Sigma
    Sigobs = set()
    for obschannel in OC:
        Sigobs = Sigobs | set(obschannel[1])
    Sigcon = auto.Sigcon
    Gtable = auto.symDict
    maxDelay,Eoci = create_CHdicts(OC)
    table = {}
    for ev in Sigma:
        if ev in Gtable:
            table.update( {ev: Gtable[ev]} )
        else:
            table.update( {ev: str(ev)} )

    _nu_ = 'nu'
    table.update( {_nu_: r'\nu'} )
    Xlist = {}
    T = []
    
    # Step 1:
    x0d = newstate( [_nu_] )
    Xd = []
    #print('--------- print(x0d); print(Xlist[x0d]); print(table[x0d]) -------------')
    #print(x0d); print(Xlist[x0d]); print(table[x0d])
    
    # Step 2:
    El = set()
    Es = set()
    for event in Sigobs:
        if event in Elo:
            label,st = rename_event(event,'_','d',Gtable,loi,ST=True)
            El = El | {label}
            table.update({label: st})
        label,st = rename_event(event,'_','r',Gtable,loi,ST=True)
        Es = Es | {label}
        table.update({label: st})
    Ei = Sigma | Es | El
    #print('------ print(El,Es,Ei) ---------')
    #print(El,Es,Ei)
    
    # Step 3:
    F = deque()
    F.append(x0d)
    
    # Step 4:
    while F:
        #print('-------------------------')
        
        # Step 4.1:
        state = F[0]
        q = Xlist[state]
        #print('------------------ q -------------------------')
        #print(q)

        # Step 4.2:
        if q == Xlist[x0d]:
            for event in Sigma:
                if event in Sigobs:
                    newStateLabel = [event]
                    newState = newstate(newStateLabel)
                    F.append(newState)
                else:
                    newState = state
                T.append( [state,event,newState] )
                #print('----------')
                #print(Xlist[state],table[event],Xlist[newState])
            #print([Xlist[x] for x in F])


        # Step 4.3:
        else:

            # Steps 4.3.1 to 4.3.5:
            isSigmaActive = True
            #Io = []; I = {}; MET = {}; rho = {}
            #iy = []; met = 0
            l = len(q)
            Il = range(0,l)
            Ilnu = [index for index in Il if q[index] in Sigobs]
            for index in Ilnu:
                if (l - index) > maxDelay[q[index]]:
                    isSigmaActive = False
            #print('isSigmaActive: ',isSigmaActive)
            #
            if isSigmaActive:
                for event in Sigma:
                    isnewstate = True
                    if event in Sigobs:
                        newStateLabel = q + [event]
                    else:
                        newStateLabel = q + [_nu_]
                        xNewList = [x for x in F if Xlist[x] == newStateLabel]
                        if xNewList:
                            isnewstate = False
                            newState = xNewList[0]
                    if isnewstate:
                        newState = newstate(newStateLabel)
                        F.append(newState)
                        #print([Xlist[x] for x in F])
                    T.append( [state,event,newState] )
                    #print('----------')
                    #print(Xlist[state],table[event],Xlist[newState])
        
            # Step 4.3.6:
            for eventset in Eoci.values():
                eventset = list(eventset)
                #print('-----------------')
                #print('eventset: ',eventset)
                Y = [i for i,j in enumerate(q) if j in eventset]
                if Y!=[]:
                    y = Y[0]
                    newStateLabel = cut(rep(q,y))
                    #print('----------')
                    #print('y,  newStateLabel: ',y,newStateLabel)
                    isnewstate = True
                    for temp in Xlist.keys():
                        if Xlist[temp] == newStateLabel:
                            newState = temp
                            isnewstate = False
                            break
                    if isnewstate:
                        newState = newstate(newStateLabel)
                        F.append(newState)
                        #print([Xlist[x] for x in F])
                    label = rename_event(q[y],'_','r',Gtable,loi)
                    T.append( [state,label,newState] )
                    #print(Xlist[state],table[label],Xlist[newState])
                    if q[y] in Elo:
                        label = rename_event(q[y],'_','d',Gtable,loi)
                        T.append( [state,label,newState] )
                        #print(Xlist[state],table[label],Xlist[newState])

        # Step 4.4:
        Xd = Xd + [state]
        #print('Xd: ',[Xlist(x) for x in Xd])
        
        # Step 4.5:
        F.popleft()
    
    # Step 5 and the construction of fsa object by using DESlab
    if loi==0:
        name = '$D$'
    else:
        name = '$D_{' + str(loi) + '}$'
    D = fsa ( Xd, Ei, T, [x0d], Xd, table = table, Sigobs = Es, Sigcon=Sigcon ,name = name)
    D.setgraphic('verifier')
    
    return D,Xlist
       
def createGaSD(G,Xs,SD,SigmaD):

    auto = G.copy()
    D,Xlist = autoD(auto,SD,SigmaD,loi=0)
    
    auto.Xm = auto.X
    D.Xm = D.X0

    mapping = []
    Er = []
    Gint = auto.copy()
    for events in auto.Sigobs:
        mapping.append((events, events + '_{r}'))
        Er.append(events + '_{r}')
    Gint = Gint.renamevents(mapping)
    
    for states in Xs:
        Gint =  Gint.deletestate(states)
    Gint = ac(Gint)
    
    Gshf = parallel(auto,D)
    Gshf.name = '$G_{shf}$'
    Gshf.setgraphic('observer')
    
    GaSD = parallel(Gshf,Gint)
    GaSD.name = '$G_a^{SD}$'
    GaSD.setgraphic('verifier')
    
    return D,Gint,Gshf,GaSD,Er

def ReleaseReach(auto,state,sigma):
    Gtemp = auto.setpar(X0 = state)
    for transition in Gtemp.transitions():
        if transition[1] not in sigma:
            Gtemp = Gtemp.deletetransition(transition)
    Gtemp = ac(Gtemp)

    return list(Gtemp.X)

def SetGamma(auto,state_list):
    event_set = []
    for states in state_list:
        event_set.extend(list(auto.Gamma(states)))

    return set(event_set)


def CSOUenfSHUFFLING(auto,D,Xs,Xu):
    
    ## busca de scc em G
    ## remove estados secretos de Gint
    ## Usa apenas linguagem e não mais SCT
    
    D.Xm = D.X0
    
    sigobs = auto.Sigobs

    flag_add_su = 0
    for states in auto.X:
        if auto.Gamma(states)== frozenset():
            if flag_add_su == 0:
                flag_add_su = 1
                auto = auto.addevent(su)
            auto = auto.addtransition([states,su,states])

    auto.Sigobs = sigobs

    xm = []
    for comp in strconncomps(auto):
        if len(comp)>1 or any(u == v for u, v in selfloop_edges(auto.Graph)):
            for states in comp:
                xm.append(states)
                
    auto = auto.setpar(Xm=xm)
    
    Gint = auto.copy()
        
    mapping = []
    Er = []
    
    for events in auto.Sigobs:
        mapping.append((events, events + '_{r}'))
        Er.append(events + '_{r}')
    Gint = Gint.renamevents(mapping)

    for states in Xs:
        Gint =  Gint.deletestate(states)
    Gint = ac(Gint)

    xm = []
    for comp in strconncomps(Gint):
        if len(comp)>1 or any(u == v for u, v in selfloop_edges(auto.Graph)):
            for states in comp:
                xm.append(states)
                
    Gint = Gint.setpar(Xm=xm)

    if flag_add_su == 1:
        auto = auto.deletevent(su)
        Gint = Gint.deletevent(su)
    
    Gshf = parallel(auto,D)
    GaSD = parallel(Gshf,Gint)
    GaSD.name = '$G_a^{SD}$'
    GaSD.setgraphic('verifier')

    Sigma_r = []
    Sigma_d = []
    for events in GaSD.Sigma:
        if '{r}' in events:
            Sigma_r.append(events)
        elif '{d}' in events:
            Sigma_d.append(events)
            
    Srd = Sigma_r + Sigma_d
    Ec = Sigma_r + Sigma_d
    GaSD.Sigcon = set(Ec)
    GaSD.Sigobs=GaSD.Sigma

##    GaSD = coac(GaSD)

    Gobs = observer(auto,list(auto.Sigobs))
    Gobs.Xm = Gobs.X

    GaSD_Obs = observer(GaSD,list(auto.Sigobs))
    GaSD_Obs.Xm = GaSD_Obs.X

    if are_langequiv(Gobs,GaSD_Obs) == True:
        print('The system is current-state opaque enforceable (CSOE)')
    else:
        print('The system is NOT current-state opaque enforceable (not CSOE)')
 
    
    V = GaSD.copy()
    V.name = '$V$'
    if Xu != []:
        
    ## Caso 1: uma vez que um estado util é alcançado, ele deve ser estimado
        StatesToBeDeleted=[]
        TransitionsToBeDeleted=[]
        for states in V.X:
            if (states[0] in Xu) and (states[2] != states[0]):        
                for transitions in V.transitions():
                    if (transitions[0] == states) and (transitions[1] in auto.Sigma):
                        TransitionsToBeDeleted.append(transitions)
##                        print('(U) causa 1:',V.symDict[transitions[0]],',',transitions[1],',',V.symDict[transitions[2]])


        ## Caso 2: Toda vez que se for estimado um estado util, a planta deve de fato estar em um estado util

            if (states[2] in Xu) and (states[0] != states[2]):         
                for transitions in V.transitions():
                    if (transitions[2] == states) and (transitions[1] in Sigma_r):
                        if states not in StatesToBeDeleted:
                            TransitionsToBeDeleted.append(transitions)
##                            print('(U) causa 2:',V.symDict[transitions[0]],',',transitions[1],',',V.symDict[transitions[2]])
                     

        for transitions in TransitionsToBeDeleted:
            V = V.deletetransition(transitions)
    V = ac(V)

    H = coac(V)
    H.name = '$H=CoAc(V)$'

    Hcopy = H.copy()


    ## (OEC 1 - OEC2)
    flag = True
    while flag == True:
        flag = False
        ## Step 3.2
        for states in H.X:
            Rreach = ReleaseReach(H,states,Srd)
            setgamma = SetGamma(H,Rreach)
            if ((auto.Gamma(states[0])& D.Gamma(states[1]))- setgamma) != frozenset():
                H = H.deletestate(states)
                flag = True
##                print('(OE) causa 2:',states[0],',',D.symDict[states[1]],',',states[2])

        H = ac(coac(H))

        
    if H == fsa():
        print (H== fsa())
        print('Utility can NOT be ensured in the system while CSO is being enforced (not UE-CSOE)')
        return GaSD,V,Hcopy,H
    else:
        O = observer(H,list(auto.Sigobs))
        O.Xm = O.X

        if are_langequiv(Gobs,O) == True:
            print('Utility can be ensured in the system while CSO is being enforced (UE-CSOE)')
        else:
            print('Utility can NOT be ensured in the system while CSO is being enforced (not UE-CSOE)')

    
    
    ## (UC 1 - UC3)
    for states in H.X:
        ## Step 4.1
        if H.Gamma(states)&set(Sigma_r) != frozenset():
            EventsToBeDeleted = H.Gamma(states)-set(Sigma_r)
            for transitions in H.transitions():
                if transitions[0]==states and (transitions[1] in EventsToBeDeleted):
                    H = H.deletetransition(transitions)
        ## Step 4.2
        if H.Gamma(states)&auto.Sigobs != frozenset():
            EventsToBeDeleted = H.Gamma(states)-auto.Sigobs
            for transitions in H.transitions():
                if transitions[0]==states and (transitions[1] in EventsToBeDeleted):
                    H = H.deletetransition(transitions)
        ## Step 4.3
        if len( H.Gamma(states)&(D.Sigma-auto.Sigma))>1:
            ## Step 4.3.1
            x_D = D.symDict[states[1]]
            i = 1
            ## Step 4.3.2
            while (( set([x_D[i-1]+'_{r}']+[x_D[i-1]+'_{d}'])& H.Gamma(states)) == frozenset()) and (i <= len(x_D)):
                i = i+1
            ## Step 4.3.3
            EventsToBeDeleted = H.Gamma(states)-set([x_D[i-1]+'_{r}']+[x_D[i-1]+'_{d}'])
            for transitions in H.transitions():
                if transitions[0]==states and (transitions[1] in EventsToBeDeleted):
                    H = H.deletetransition(transitions)
    H = ac(H)
    H.name = '$H_{ue}$'

    O = observer(H,list(auto.Sigobs))
    O.Xm = O.X
    
    if are_langequiv(Gobs,O) == True:
        print('P(L(H_{ue})) = L(G)')
    else:
        print('P(L(H_{ue})) != L(G)')

    
    return GaSD,V,Hcopy,H

def cso_shuffle_deletion_function(G, SD, SigmaD, Xs, Xu):
    """
    This function returns an automaton that shuffle and/or delete events to 
    enforce current state opacity to an automaton G, where SD is the maximum time 
    a event can be hold until it needs to be released, SigmaD the events that can be deleted,
    Xs are the secret states and Xu are .
    If it returns an empty automaton means the enfocement can be done.
    
    ---------------
    X = [ 0, 1, 2, 3,4,5,6,7,8 ]
    E = [ a, b, c]
    T = [ (0,a,1), (0,b,6), (1,b,2), (1,c,4), (2,c,3), (4,b,5), (6,c,7), (7,a,8)]
    X0 = [ 0 ]
    Xs = [ 3 ]
    Xu = [ 4 ]
    G = fsa ( X , E , T , X0, name = '\$G\$')
    SD = [(2,[a]), (0,[b]), (1,[c])]
    SigmaD = [ ]
    
    sdf = cso_shuffle_deletion_function(G, SD, SigmaD, Xs, Xu)
    draw(G, sdf, 'figure')
    """
    D,Gint,Gshf,GaSD,Er = createGaSD(G,Xs,SD,SigmaD)

    GaSD_new,V,H,Huc = CSOUenfSHUFFLING(G,D,Xs,Xu)
    return Huc

"""
syms('0 1 2 3 4 5 6 7 8 9 10 11 a b c d su')

table = [(0,'0'),(1,'1'),(2,'2'),(3,'3'),(4,'4'),(5,'5'),(6,'6'),(7,'7'),(8,'8'),(9,'9'),(10,'10'),(11,'11'),(a,'a'),(b,'b'),(c,'c'),(d,'d'),(su,r'\sigma_u')]
  
X2 = [ 0, 1, 2, 3,4,5,6,7,8 ] # state set
E2 = [ a, b, c] # event set
T2 = [ (0,a,1), (0,b,6), (1,b,2), (1,c,4), (2,c,3), (4,b,5), (6,c,7), (7,a,8)] # transitions
X02 = [ 0 ] # set of initial states
Xs2 = [ 3 ] # set of secret states
Xu2 = [ 4 ] # set of utility states
Xm2 = Xs2 # set of marked states
##Eo = [a, b, c]
G2 = fsa ( X2 , E2 , T2 , X02 , [ ] , table = table , name = '$G_2$') #Building G
SD2 = [(2,[a]), (0,[b]), (1,[c])] #set if Step Delay bounds
SigmaD2 = [ ] #deletable events set


D2,Gint2,Gshf2,GaSD2,Er2 = createGaSD(G2,Xs2,SD2,SigmaD2)
print('------------------------------------')

GaSD_new,V2,H2,Huc2 = CSOUenfSHUFFLING(G2,D2,Xs2,Xu2)
             
"""
  