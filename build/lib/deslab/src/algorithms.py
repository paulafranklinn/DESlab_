"""Mathematical functions for the automata objects"""
#    Copyright (C) 2011-2012 by
#    Leonardo Bermeo <lbermeoc@unal.edu.co>
#    Joao Carlos Basilio <basilio@poli.ufrj.br>
#    GNU license.
# Cj3Uk7yV5JS8


import networkx as nx
from deslab.src.automatadefs import *
from deslab.src.def_const import *
import copy
from deslab.src.structure import *
from deslab.graphics.drawing import graphic

class compCount:
    """ This class is a counter for dump states created in some
    operations (for instance, complement). When we instantiate
    a new automaton and call a new counter, that counter is 
    incremented by 1.
    """
    counter = -1 # static variable, counts no. of instances
    def __init__(self):        
        compCount.counter += 1


    
def dfs(H, source=None):
    """Produce edges in a depth-first-search starting at source.
    it is the base of many algorithms used in accesibility"""
    from deslab.src.automatadefs import fsa
    if isinstance(H,fsa):
        G = H.Graph
    elif isinstance(H, nx.MultiDiGraph):
        G = H
    nodes=list(source)
    visited=set()
    for start in nodes:
        if start in visited:
            continue
        visited |= set([start])
        stack = [(start,iter(G[start]))]
        while stack:
            parent,children = stack[-1]
            try:
                child = next(children)
                if child not in visited:                    
                    visited |= set([child])                    
                    stack.append((child,iter(G[child])))
            except StopIteration:
                stack.pop()
    visited = frozenset(visited)
    return visited



def ac(automaton):
    """ This function calculates the accessible part
    of an automaton object. """     
    from deslab.src.automatadefs import fsa, create_FSA_transdicts
    auto = automaton.copy() 
    X0,X,Sigma,Xm = auto.X0,auto.X,auto.Sigma,auto.Xm     
    if (auto.empty==True)|(X0 == EMPTYSET):        
        return fsa()    
    Xac = dfs(auto.Graph,X0) #calculating accesibility  
    Xmac = Xm & Xac    
    if Xac == set([]):
        return fsa()
    auto.X = Xac 
    auto.Graph.remove_nodes_from(X-Xac)  
    auto.gammaDict, auto.deltaDict, auto.infoDict = create_FSA_transdicts(auto.Graph,Xac,Sigma,X0)
    auto.Xm = Xmac
    return auto


def coac(automaton):  
    """ This function calculates the coaccessible part
    of an automaton object. """  
    from deslab.src.automatadefs import create_FSA_transdicts,fsa
    auto=automaton.copy()
    X0,X,Sigma,Xm = auto.X0,auto.X,auto.Sigma,auto.Xm    
    if (auto.empty==True)  |  ( Xm == EMPTYSET):
        return fsa()    
    GraphReversed = auto.Graph.reverse()
    Xcoac = set([])    
    Xcoac = dfs(GraphReversed,Xm)  
    auto.Graph.remove_nodes_from(X-Xcoac)  # we left coaccesible part
    auto.gammaDict, auto.deltaDict, auto.infoDict = create_FSA_transdicts(auto.Graph,Xcoac,Sigma,X0)    
    if (X0 & Xcoac) == EMPTYSET: 
        # if initial state isn't in Xcoac result is empty
        return fsa()    
    auto.X = Xcoac
    auto.Xm = Xm
    auto.X0 = X0 & Xcoac
    return auto

def trim(automaton):
    """ This function calculates a trimmed 
    version of an automaton object G passed to it. """
    trimmed=coac(automaton)
    trimmed=ac(trimmed)
    trimmed.setgraphic('normal')
    return trimmed

def invproj(self,sigma):
    """ This function calculates  the inverse projection 
    of the automaton object G  """
    from deslab.src.automatadefs import create_FSA_transdicts
    
    if not self.Sigma <= set(sigma):
        raise invalidArgument("The list of events 'sigma' must contain all of the events of the input automaton")

    auto = self.copy()
    X , Gamma, Sigma  = auto.X, auto.Gamma, auto.Sigma
    X0, name = auto.X0, auto.name
    sigmaEuo = set(sigma) - Sigma
    # case the observed set is the same    
    if (sigmaEuo == []) | (sigmaEuo == EMPTYSET) :
        return self
    for state in X:
        for event in sigmaEuo:
            auto.Graph.add_edge(state,state,key=event,label=event) # adding transitions
    sigmas = Sigma | set(sigma)
    auto.gammaDict, auto.deltaDict, auto.infoDict=create_FSA_transdicts(auto.Graph, X, sigmas, X0)     
    auto.name = 'untitled'
    auto.Sigma = Sigma | set(sigma)
    auto.Sigobs = Sigma
    return auto


def pclosure(self):
    """This function receives an automaton as argument and returns
    an automaton which marks the prefix of the language marked
    by the received automaton."""
    auto = self.copy() # copy of automaton
    auto = trim(auto)
    auto = auto.setpar(Xm = auto.X)       
    return auto

def langquotient(H1,H2):
    """
    This function computes automaton M such that Lm(M) = Lm(H1) / Lm(H2),
    where / denotes the (right) quotient of Lm(H1) and Lm(H2)
    """
    
    if len(H1.Xm)==0 or len(H2.Xm)==0:
        return fsa()
    elif len(H2.X)==1:
        H1x = H1.copy()
        SetOfKleenClosure = H2.Gamma(list(H2.X)[0])
        Sigma_rem = H1x.Sigma - SetOfKleenClosure
        for x in H1x.X:
            for event in H1x.Gamma(x) & Sigma_rem:
                H1x =  H1x.deletetransition([x,event,H1x.delta(x,event)])
        H1x = H1x.setpar(X0 = H1x.X)
        H1x = coac(H1x)
        Xm = H1x.X
    else:
        Xm = []
        for x in H1.X:
            H1x = H1.copy()
            H1x = H1x.setpar(X0 = [x])
            auto = H1x & H2
            if auto.Xm:
                Xm = Xm + [x]
    M = H1.copy()
    M = M.setpar(Xm=Xm)
    M = trim(M)
    M = M.setpar(Sigobs= H1.Sigobs)
    M = M.setpar(Sigcon= H1.Sigcon)
    M = M.setpar(name= 'quotient(' + H1.name + ',' + H2.name + ')')
    return M

def langdiff(self,other):
    """ This function calculates the set difference between
    the languages marked by the automata G1 and G2.
    It returns the automaton D such that Lm(D) = Lm(G2)\Lm(G1)
    
    Example
    -------
    
    syms('q1 q2 q3 a1 b1 e f')
    simtab=[(q1,'q_1'),(q2,'q_2'),(q3,'q_3'),
            (a1,'a_1'),(b1,'b_1')]
    X = [q1,q2,q3]
    Sigma = [a1,b1,e]
    X0 = [q1]
    Xm = [q2]
    Xm2 = [q2,q3]
    Sigma2= [a1,b1, f]
    T =[(q1,a1,q2),(q2,b1,q3)]
    T2 = [(q1,a1,q2),(q2,b1,q3)]
    G1 = fsa(X,Sigma,T,X0,Xm,simtab,name='$G_1$')
    G2 = fsa(X,Sigma,T2,X0,Xm2,simtab,name='$G_2$')
    D=langdiff(G2,G1)
    draw(G1,G2,D)
    """   
    # imports from modules  
    from deslab.src.automatadefs import fsa
    
    if self == fsa(): # case first automaton is empty
        return fsa()
    elif other == fsa(): # case 2nd automaton is empty
        return self   
        
    lang_diff = self & complement(other)
    lang_diff = trim(lang_diff)
    #lang_diff.setpar(name='untitled')#,type='plant')
    return lang_diff

def complement(self):       
    """ This function calculates the complement of 
    the language marked by the automaton G. It returns
    the automaton C such that Lm(C) = Lm(G)^C
    
    Example
    -------
    syms('q1 q2 q3 a1 b1 e f')
    table = [(a1,'a'),(b1,'b_1'),(q1,'q_1'),(q2,'q_2'),(q3,'q_3')]
    X = [q1,q2,q3]
    Sigma = [a1,b1,e]
    X0 = [q1]
    Xm = [q2]
    T =[(q1,a1,q2),(q2,b1,q3)]
    G1 = fsa(X,Sigma,T,X0,Xm,table,name='$G_1$ -- example 1')
    Gc=complement(G1)
    draw(G1,Gc)
    """            
    from deslab.src.automatadefs import create_FSA_transdicts    
    # we create dump states with different labels 
    auto = self.copy()
    # we create a unique label for dump_state
    countInstance=compCount()
    count=countInstance.counter  
    dump = '_XD_'+str(count) # dump state
    dump_latex = 'X_D'# latex string of dump state 
    
    X , Gamma, Sigma  = auto.X, auto.Gamma, auto.Sigma
    X0, name = auto.X0, auto.name
     
    auto.Graph.add_node(dump,label='s'+str(auto.Graph.order()+1))
    auto.X = X | frozenset([dump]) # new state set        
 
    for event in Sigma:
        auto.Graph.add_edge(dump,dump,key=event,label=event)  
    
    for state in X:   
        validSet = Sigma - Gamma(state) 
        for event in validSet: 
            auto.Graph.add_edge(state,dump,key=event,label=event)
    auto.gammaDict, auto.deltaDict, auto.infoDict = create_FSA_transdicts(auto.Graph,X,Sigma,X0)
    auto = auto.setpar(Xm = auto.X - auto.Xm)
    auto.symDict[dump] = dump_latex 
    auto = trim(auto)  
    return auto


def complete(self):  
     
    """ This function calculates the complete automaton
    that generates Sigma^*, where Sigma is the input 
    alphabet. The output automaton marks the same language
    as  automaton G passed as input argument.
    
    """            
    from deslab.src.automatadefs import create_FSA_transdicts    
    # we create dump states with different labels 
    auto = self.copy()
    X , Gamma, Sigma  = auto.X, auto.Gamma, auto.Sigma
    X0, name = auto.X0, auto.name
    dump = 'D' # dump state
    dump_latex = 'X_D'# latex string of dump state
      
    auto.Graph.add_node(dump,label='s'+str(auto.Graph.order()+1))
    auto.X = X | frozenset([dump]) # new state set        
 
    for event in Sigma:
        auto.Graph.add_edge(dump,dump,key=event,label=event)  
    
    for state in X:   
        validSet = Sigma - Gamma(state) 
        for event in validSet: 
            auto.Graph.add_edge(state,dump,key=event,label=event)
    auto.gammaDict, auto.deltaDict, auto.infoDict = create_FSA_transdicts(auto.Graph,X,Sigma,X0)
    auto.name = 'untitled'   
    auto.symDict[dump] = dump_latex 
    return auto


def proj(self, Sigma_o = UNDEFINED):
    Proj = observer(self,Sigma_o)
    Proj.graphic = graphic('normal')
    Proj = Proj.renamestates('number')
    return Proj
    
def sigmakleeneclos(sigma,x0 ='s0',label = 's_0'):
    
    """ This function calculates a single state automaton
    that marks and generates the language Sigma^*, where
    Sigma is the set (list) of events passed as argument.
    The name of the single state can be set with the argument
    x0.
    """    
    from deslab.src.automatadefs import fsa    
    if isinstance(sigma,list):
        Sigma = frozenset(sigma)
    elif isinstance(sigma,set):
        Sigma = frozenset(sigma)
    elif isinstance(sigma,frozenset):
        Sigma = sigma
    else :
        raise invalidArgument('Sigma must be defined using set, frozenset or list')
    
    X0 = frozenset([x0])
    X, Xm = X0, X0
    transition = []
    for event in Sigma:
        transition.append([x0,event,x0])        
    label = [(x0,label)]
    auto = fsa(X,Sigma,transition,X0,Xm,table=label) 
    return auto

def union(self, other, dfa = True):
    """
    union(G1,G2,dfa)
    This function calculates an automaton that marks
    the language Lm(G1)+Lm(G2), where G1 and G2 are the input
    automata. The resulting automaton is a deterministic automaton
    that marks the union of the languages marked by the input
    automata.
    
    Example
    -------
    
    syms('s1 s2 s3 s4 a b c d e f')
    X1=[s1,s2,s3,s4] #definindo estados
    E1=[a,b,c] #definindo eventoson
    Transition1=[(s1,a,s2),(s2,b,s3),(s3,c,s4)]
    X0=[s1]
    Xm=[s2,s3,s4]
    G1=fsa(X1,E1,Transition1,X0,Xm,name='$G_1$')
    X2=[s1,s2,s3,s4] #definindo estados
    E2=[d,e,f] #definindo eventoson
    Transition2=[(s1,d,s2),(s2,e,s3),(s3,f,s4),(s1,e,s3),(s1,f,s4),(s2,f,s4)]
    X0,Xm=[s1],[s1,s2]
    G2=fsa(X2,E2,Transition2,X0,Xm, name='$G_2$')
    U=union(H1,H2)
    draw(G1,G2,U)
    
    """  
    # imports from modules  
    from deslab.src.automatadefs import fsa
    from deslab.src.comparison import isitempty,isitemptymarked  

    # G1 marks an empty language and G2 does not
    if isitemptymarked(self) and not isitemptymarked(other):
        return other
    # G1 and G2 marks empty languages
    elif isitemptymarked(self) and isitemptymarked(other):
        return fsa()
    # G1 marks a nonempty language and G2 do.
    elif not isitemptymarked(self) and isitemptymarked(other):
        return self

    unionaut=~(trim((~self)&(~other)))

    tex=[]
    tex.append((list(unionaut.X0)[0],0))
    states = list(unionaut.X-unionaut.X0)

    for i in range(len(states)):
        tex.append((states[i],i+1))

    unionaut=unionaut.renamestates(tex)
    
    return unionaut 


def concatenation(self, other, dfa = True):
    """ This function calculates an automaton that marks
    the language Lm(G1)Lm(G2), where G1 and G2 are the input
    automata. The resulting automaton is a deterministic automaton
    that marks the concatenation of the languages marked by the input
    automata. 
    
    Example
    -------
    syms('q1 q2 q3 a1 b1 e f')
    table = [(a1,'a_1'),(b1,'b_1'),(q1,'q_1'),(q2,'q_2'),(q3,'q_3')]
    X = [q1,q2,q3]
    Sigma = [a1,b1,e]
    X0 = [q1]
    Xm = [q1,q2]
    T =[(q1,a1,q2),(q2,b1,q3)]
    G1 = fsa(X,Sigma,T,X0,Xm,table,name='$G_1$')
    Xm2 = [q2,q3]
    Sigma2= [a1,b1, f]
    T2 = [(q1,a1,q2),(q2,b1,q3)]
    G2 = fsa(X,Sigma,T2,X0,Xm2,table,name='$G_2$')
    C= concatenation(G1,G2)
    draw(G1,G2,C)
    """     
    from deslab.src.automatadefs import fsa    
    # inner functions used in union operation   
    def rename(state_set,name):
        """This is a simple function for renaming the states
        of automaton to minimize name conflicts"""
        renamed = [(state,name) for state in state_set]
        return renamed
    def generate_table():
        """this function updates the symbolic dictionary
        of the automata resulting of concatenation in order
        to preserve the labels of events of the input
        automata"""
        table = {}
        for event in Sigma_c:
            if event in self.symDict:
                table.update({event:self.symDict[event]})
            if event in other.symDict:
                table.update({event:self.symDict[event]})
        return table
    # starting main program  
    # renaming with easy names 

    if self == fsa():
        return fsa()
    elif other == fsa():
        return fsa()        
    
    Graph1, Graph2 = self.Graph, other.Graph
    X01, X02 = self.X0, other.X0
    Sigma1, Sigma2 = self.Sigma, other.Sigma
    X1, X2 = self.X, other.X
    Xm1, Xm2 = self.Xm, other.Xm  
    # initializing parameters
    transition=[]
    countInstance2 = compCount()    
    new_name1 = countInstance2.counter    
    countInstance2 = compCount()
    new_name2 = countInstance2.counter 
    #adding transitions of G1
    for edge in Graph1.edges(keys=True):     
        x_a = (edge[0], new_name1)
        x_b = (edge[1], new_name1)        
        transition.append([x_a, edge[2], x_b])
    # adding transition of G2
    for edge in Graph2.edges(keys=True):     
        x_a = (edge[0], new_name2)
        x_b = (edge[1], new_name2)        
        transition.append([x_a, edge[2], x_b])
    # connecting final states of G1 with initial states of G2
    for xm_1 in Xm1 :
        for x0_2 in X02 :
            transition.append([(xm_1,new_name1), epsilon, (x0_2, new_name2)])
    # renaming initial states
    X0_cren = rename(X01,new_name1) 
    X0_c = frozenset(X0_cren)
    # renaming states     
    X1_ren, X2_ren = rename(X1,new_name1), rename(X2,new_name2)
    X_c = X1_ren + X2_ren 
    # renaming marked states
    Xm_c =  rename(Xm2,new_name2)
    # building alphabet of the concatenation automaton
    Sigma_c = Sigma1  | Sigma2 | EMPTYSTRINGSET 
    table_c =  generate_table()       
    concat = fsa(X_c, Sigma_c, transition, X0_c, Xm_c,
                 table=table_c, name='untitled')    
    Sigobs_c = self.Sigobs | other.Sigobs 
    Sigcon_c = self.Sigcon | other.Sigcon
    concat = concat.renamestates('number') 
    if dfa :
        concat = epsilonobserver(concat)
        concat = concat.renamestates('number')  
        concat = concat.setpar(Sigcon=Sigcon_c, Sigobs=Sigobs_c)
        concat.setgraphic('normal', direction='LR')         
    return concat 
     
     



def observer(self, Sigma_o=UNDEFINED):    
    """ This function calculates a deterministic automaton
    which marks and generates the projection of the marked and
    generated languages of the input automaton G over the alphabet
    of observable events Sigma_o, defined by the user. If a list (or set)
    of observable events is not provided, the function observer takes
    the set Sigobs of the input automaton as the set of observable
    events"""
    from deslab.src.automatadefs import fsa
    # inner functions    
    def Gamma_obs(set): # Gamma observer calcularion                        
        active_events = EMPTYSET
        for q in set:
            active_events = active_events | self.Gamma(q)
        active_events = active_events & Sigma_obs
        return active_events  

    def latexname(Xnew): 
        """This function takes a frozenset of the observer 
        automaton and generate a latex label, based
        on the dictionaries of the input automaton""" 
        Xnew = list(Xnew)# in observer objects
        #print(Xnew)
        #Xnew.sort()
        name_state=''            
        for xn in Xnew:
            if xn in self.symDict:
                name_state = name_state+self.symDict[xn]+',' 
            else:
                name_state = name_state+str(xn)+','       
        name_state='\{'+name_state.rstrip(',')+'\}'          
        return name_state
    # Main program    
    # easy name for parameters
    X0, X = self.X0, self.X
    Sigma, Xm, Sigobs = self.Sigma, self.Xm, self.Sigobs    
    # determinig the set of  unobservable events
    if (Sigma_o == UNDEFINED) : # Sigma_o was not defined
        Sigma_obs = Sigobs
        Sigma_uobs = Sigma - Sigobs        
    elif isinstance(Sigma_o,list) | isinstance(Sigma_o,set):      
        Sigma_obs = frozenset(Sigma_o)
        Sigma_uobs = Sigma - Sigma_obs         
    elif isinstance(Sigma_o,frozenset):
        Sigma_obs = Sigma_o 
        Sigma_uobs = Sigma - Sigma_o
    else: 
        raise invalidArgument('Observed set must be of type set, frozenset or list')
    # in the case of not observable events
    # we return the input automaton         
    if (self.is_dfa() &  (Sigma_obs == self.Sigma)):
        return self       
    if self.infoDict['hasEpsilon'] and Sigma_obs != self.Sigma - EMPTYSTRINGSET:
        raise invalidArgument('Try epsilonobserver(G) instead')
    #initial values of parameters
    X0_obs = self.unobsreach(X0,Sigma_obs)    
    X_obs = [X0_obs]
    S = [X0_obs] # unobs. reach of X0 as initial value of stack
    transition = []
    Xm_obs = []   
    # updating X0_p label
    table_obs = {X0_obs: latexname(X0_obs)}
    # testing of X0obs is marked
    if X0_obs & Xm != EMPTYSET: 
        Xm_obs.append(X0_obs)
    # main stack cycle      
    while S != []: # beginning of subset construction      
        Q = S.pop()                     
        for sigma in Gamma_obs(Q):
            Q_next = self.deltaobs(Q,sigma,Sigma_obs) # reached set 
            if Q_next not in X_obs: # net state found
                X_obs.append(Q_next)
                S.append(Q_next) 
                if Q_next & Xm != EMPTYSET: # marking that state
                    Xm_obs.append(Q_next) 
                table_obs.update({Q_next: latexname(Q_next)})  # creating the latex table
            transition.append([Q,sigma,Q_next])
            if sigma in self.symDict:
                table_obs.update({sigma:self.symDict[sigma]}) 
    # defining the observer DFA           
    observer = fsa(X_obs, Sigma_obs, transition, [X0_obs],
                   Xm_obs, Sigobs = Sigma_obs, table=table_obs, Sigcon=self.Sigcon) 
    observer.graphic = graphic('observer') # graphic property for square states
    return observer 


def parallelnondet(self,other,simplify):
    """This function builds the automaton corresponding
    to the parallel composition of the input automata.
    Specifically, parallelnondet assumes that the input
    automata (one or both of them) are non deterministic """
    from deslab.src.automatadefs import fsa
    
    def pars(var):
        resp=()
        for i in var:
            if type(i)==tuple:
                resp+=pars(i)
            else:
                resp+=(i,)
        return resp

    def Gamma_p(p):
        """definition of the Gamma function of parallel
        composition of the input automata for the tuple
        p = (x1,x2)"""
        x1, x2 = p # defined for the tuple p
        # these are the active events
        active_events = (Gamma1(x1) & Gamma2(x2)) | (Gamma1(x1) - Sigma2)
        active_events = active_events | (Gamma2(x2) - Sigma1) 
        return active_events
    
    def delta_p(p,sigma):
        """definition of the delta transition  function
         of parallel composition of the input automata. It is
         defined for the tuple p = (x1,x2) and event sigma
         the output is a set of states that can be reached
         from p"""
        x1,x2 = p   
                      
        if sigma in Gamma1(x1) & Gamma2(x2) :
            X1_n, X2_n = delta1(x1,sigma),delta2(x2,sigma)                     
        elif sigma in Gamma1(x1)-Sigma2 :
            X1_n, X2_n = delta1(x1,sigma), frozenset([x2])             
        elif sigma in Gamma2(x2)-Sigma1 :
            X1_n, X2_n = frozenset([x1]), delta2(x2,sigma)          
        else :
            X1_n, X2_n = EMPTYSET, EMPTYSET                       
        # it returns a tuple of sets because the function 
        # parallelnonder considers input automata as non
        # deterministic
        return X1_n, X2_n 

    
    def latexname(p,simplify):  
        """This function takes a tuple of the composite 
        automaton and generate a latex pretty name, based
        on the dictionaries of each input automaton"""
        x1, x2 = p
        if x1 in self.symDict:
            x1_name = str(self.symDict[x1])
        else: 
            x1_name = str(x1)
        if x2 in other.symDict:
            x2_name = str(other.symDict[x2])
        else: 
            x2_name = str(x2)
        # pretty tuple print as, for instance, (\pi,\Gamma)
        name = '('+x1_name+','+x2_name+')'
        if simplify:
            name = name.replace("(","")
            name = name.replace(")","")
            name = '('+name+')'
        return name                           
    # Main program 
    # easy names for inner variables
    X1, X2 = self.X, other.X
    X01, X02 = self.X0, other.X0
    Xm1, Xm2 = self.Xm, other.Xm 
    Sigma1, Sigma2 = self.Sigma, other.Sigma
    Gamma1, Gamma2 = self.Gamma, other.Gamma 
    delta1, delta2 = self.__delta__, other.__delta__    
    Sigobs1,Sigobs2 = self.Sigobs, other.Sigobs         
    Sigcon1,Sigcon2 = self.Sigcon, other.Sigcon  
    
    # initial values for the composite automaton
    X0_p = [(x01,x02) for x01 in X01 for x02 in X02]
    S = [(x01,x02) for x01 in X01 for x02 in X02]
    X_p = frozenset(X0_p)
    Xm_p = EMPTYSET
    Sigobs_p = Sigobs1|Sigobs2 
    Sigcon_p = Sigcon1|Sigcon2
    Sigma_p = Sigma1 | Sigma2
    transition=[]
    # generating labeling table for initial states
    table_p = {}
    for p in X0_p: # labeling initial states
        table_p.update({p: latexname(p,simplify)})
        p1,p2 = p
        if (p1 in Xm1) & (p2 in Xm2):
            Xm_p = Xm_p | frozenset([p])
        
    # main stack cycle     
    while S != [] :
        p = S.pop(); # taking a tuple     
        # evaluating in tha active event set
        for sigma in Gamma_p(p): 
            Q1_n, Q2_n = delta_p(p, sigma)
            Q = frozenset((q1,q2) for q1 in Q1_n  for q2 in Q2_n)
            # Q contains the tuples fo reached states
            for q in Q:                 
                if q not in X_p: 
                    # updating states of composite automaton and stack                                                    
                    X_p = X_p | frozenset([q])           
                    S.append(q)
                    q1, q2 = q
                    # updating marked states
                    if (q1 in Xm1) & (q2 in Xm2):
                        Xm_p =Xm_p | frozenset([q]) 
                    #updating latex dictionary
                    table_p.update({q: latexname(q,simplify)})
                transition.append([p,sigma,q])  
            if sigma in self.symDict:
                table_p.update({sigma: self.symDict[sigma]}) 
            elif sigma in other.symDict:
                table_p.update({sigma: other.symDict[sigma]})  
    #Simplify
    if simplify:
        X_p2=()
        for i in X_p:
            X_p2+=(pars(i),)
        X_p = X_p2
        for i in range(len(transition)):
            for j in range(len(transition[i])):
                if type(transition[i][j])==tuple:
                    transition[i][j]=pars(transition[i][j])
        for i in range(len(X0_p)):
            if type(X0_p[i])==tuple:
                X0_p[i]=pars(X0_p[i])
        Xm_p=list(Xm_p)
        for i in range(len(Xm_p)):
            if type(Xm_p[i])==tuple:
                Xm_p[i]=pars(Xm_p[i])
        Xm_p=frozenset(Xm_p)
        keys=list(table_p.keys())
        values=list(table_p.values())
        for i in range(len(keys)):
            if type(keys[i])==tuple:
                keys[i]=pars(keys[i])
        table_p={}
        for i in range(len(values)):
            if values[i].count(")")>0:
                values[i]="".join(values[i].split("("))
                values[i]="".join(values[i].split(")"))
                values[i]="("+values[i]+")"
            table_p[keys[i]]=values[i]

    # finally, we build the composite automaton                                
    parnondet = fsa(X_p, Sigma_p, transition, X0_p, Xm_p, table=table_p,
                  Sigobs = Sigobs_p, Sigcon=Sigcon_p,name='untitled')    
    return parnondet



def paralleldet(self,other,simplify):
    """This function builds the automaton corresponding
    to the parallel composition of the input automata.
    Specifically, paralleldet is optimized for deterministic
    automata and it produces an error if the user tries to
    use it with nondeterministic automata"""
    from deslab.src.automatadefs import fsa

    def pars(var):
        resp=()
        for i in var:
            if type(i)==tuple:
                resp+=pars(i)
            else:
                resp+=(i,)
        return resp

    def Gamma_p(p):
        """ definition of the Gamma function of parallel
        composition of the input automata for the tuple
        p = (x1,x2)"""
        x1, x2 = p
        active_events = (Gamma1(x1) & Gamma2(x2)) | (Gamma1(x1) - Sigma2)
        active_events = active_events | (Gamma2(x2) - Sigma1)   
        return active_events
    
    def delta_p(p,sigma):
        """ definition of the delta transition  function
        of parallel composition of the input automata. It is
        defined for the tuple p = (x1,x2) and event sigma
        the output is the single state that can be reached
        from p
        """
        x1,x2 = p                 
        if sigma in Gamma1(x1) & Gamma2(x2) :
            x1_n, x2_n = delta1(x1,sigma),delta2(x2,sigma)                     
        elif sigma in Gamma1(x1)-Sigma2 :
            x1_n, x2_n = delta1(x1,sigma), x2             
        elif sigma in Gamma2(x2)-Sigma1 :
            x1_n, x2_n = x1, delta2(x2,sigma)          
        else :
            x1_n, x2_n = EMPTYSET, EMPTYSET
        return x1_n, x2_n 

    
    def latexname(p,simplify):  
        """This function takes a tuple of the composite 
        automaton and generate a latex pretty name, based
        on the dictionaries of each input automaton"""
        x1, x2 = p
        if x1 in self.symDict:
            x1_name = self.symDict[x1]
        else: 
            x1_name = str(x1)
        if x2 in other.symDict:
            x2_name = other.symDict[x2]
        else: 
            x2_name = str(x2)

        name = '('+x1_name+', '+x2_name+')'
        if simplify:
            name = name.replace("(","")
            name = name.replace(")","")
            name = '('+name+')'
        return name    
                            
        

        return new_name
    # Main program
    # easy names for inner variables
    X1, X2 = self.X, other.X
    X01, X02 = self.X0, other.X0
    Xm1, Xm2 = self.Xm, other.Xm 
    Sigma1, Sigma2 = self.Sigma, other.Sigma
    Gamma1, Gamma2 = self.Gamma, other.Gamma 
    delta1, delta2 = self.delta, other.delta   
    Sigobs1,Sigobs2 = self.Sigobs, other.Sigobs         
    Sigcon1,Sigcon2 = self.Sigcon, other.Sigcon  
    # initial values for composite automaton
    p0 = (list(X01)[0], list(X02)[0])
    X0_p = [p0]
    S = [p0]
    X_p = frozenset(X0_p)
    Xm_p = EMPTYSET
    Sigobs_p = Sigobs1|Sigobs2 
    Sigcon_p = Sigcon1|Sigcon2
    Sigma_p = Sigma1 | Sigma2
    transition=[]
    table_p = {}
    table_p.update({p0: latexname(p0,simplify)})
    p01, p02 = p0
    # is initial state marked ?
    if (p01 in Xm1) & (p02 in Xm2):                
        Xm_p = Xm_p | frozenset([p0])     
    # main stack cycle    
    while S != [] :
        p = S.pop(); # taking a tuple     
        for sigma in Gamma_p(p):
            q = delta_p(p, sigma) # reached state                                 
            if q not in X_p:
                # updating state and stack                                                     
                X_p = X_p | frozenset([q])           
                S.append(q)
                q1, q2 = q
                if (q1 in Xm1) & (q2 in Xm2):
                    # updating marked states
                    Xm_p =Xm_p | frozenset([q])
                # updating labeling table 
                table_p.update({q: latexname(q,simplify)})                
            transition.append([p,sigma,q])  
            if sigma in self.symDict:
                table_p.update({sigma: self.symDict[sigma]}) 
            elif sigma in other.symDict:
                table_p.update({sigma: other.symDict[sigma]})  
    # building the composite automaton
    if simplify:
        X_p2=()
        for i in X_p:
            X_p2+=(pars(i),)
        X_p = X_p2
        for i in range(len(transition)):
            for j in range(len(transition[i])):
                if type(transition[i][j])==tuple:
                    transition[i][j]=pars(transition[i][j])
        for i in range(len(X0_p)):
            if type(X0_p[i])==tuple:
                X0_p[i]=pars(X0_p[i])
        Xm_p=list(Xm_p)
        for i in range(len(Xm_p)):
            if type(Xm_p[i])==tuple:
                Xm_p[i]=pars(Xm_p[i])
        Xm_p=frozenset(Xm_p)
        keys=list(table_p.keys())
        values=list(table_p.values())
        for i in range(len(keys)):
            if type(keys[i])==tuple:
                keys[i]=pars(keys[i])
        table_p={}
        for i in range(len(values)):
            if values[i].count(")")>0:
                values[i]="".join(values[i].split("("))
                values[i]="".join(values[i].split(")"))
                values[i]="("+values[i]+")"
            table_p[keys[i]]=values[i]
    pardet = fsa(X_p, Sigma_p, transition, X0_p, Xm_p, table=table_p, Sigobs = Sigobs_p, Sigcon=Sigcon_p,name='untitled')
    return pardet



def productnondet(self,other,simplify):
    """This function builds the automaton corresponding
    to the product composition (intersection) of the input automata.
    Specifically, productnondet assumes that the input
    automata (one or both of them) are not deterministic """
    from deslab.src.automatadefs import fsa

    def pars(var):
        resp=()
        for i in var:
            if type(i)==tuple:
                resp+=pars(i)
            else:
                resp+=(i,)
        return resp
    
    def Gamma_prod(p):
        """definition of the Gamma function of product
        composition of the input automata for the tuple
        p = (x1,x2)"""
        x1, x2 = p # defined for the tuple p
        # these are the active events
        active_events = Gamma1(x1) & Gamma2(x2)
        return active_events
        
    def latexname(p,simplify):  
        """This function takes a tuple of the composite 
        automaton and generate a latex pretty name, based
        on the dictionaries of each input automaton"""
        x1, x2 = p
        if x1 in self.symDict:
            x1_name = self.symDict[x1]
        else: 
            x1_name = x1
        if x2 in other.symDict:
            x2_name = other.symDict[x2]
        else: 
            x2_name = x2
        # pretty tuple print as, for instance, (\pi,\Gamma)
        name = '('+str(x1_name)+','+str(x2_name)+')'
        if simplify:
            name = name.replace("(","")
            name = name.replace(")","")
            name = '('+name+')'
        return name                           
    # Main program 
    # easy names for inner variables
    from deslab.src.comparison import isitempty
    from deslab.src.automatadefs import fsa
 
    # G1 marks an empty language and G2 does not
    if self.empty and not other.empty:
        return fsa()
    # G1 and G2 marks empty languages
    elif self.empty and other.empty:
        return fsa()
    # G1 marks a nonempty language and G2 do.
    elif not self.empty and other.empty:
        return fsa() 
    X1, X2 = self.X, other.X
    X01, X02 = self.X0, other.X0
    Xm1, Xm2 = self.Xm, other.Xm 
    Sigma1, Sigma2 = self.Sigma, other.Sigma
    Gamma1, Gamma2 = self.Gamma, other.Gamma 
    delta1, delta2 = self.__delta__, other.__delta__    
    Sigobs1,Sigobs2 = self.Sigobs, other.Sigobs         
    Sigcon1,Sigcon2 = self.Sigcon, other.Sigcon      
    # initial values for the composite automaton
    X0_p = [(x01,x02) for x01 in X01 for x02 in X02]
    S = [(x01,x02) for x01 in X01 for x02 in X02]
    X_p = frozenset(X0_p)
    Xm_p = EMPTYSET
    Sigobs_p = Sigobs1|Sigobs2 
    Sigcon_p = Sigcon1|Sigcon2
    Sigma_p = Sigma1 | Sigma2
    transition=[]
    # generating labeling table for initial states
    table_p = {}
    for p in X0_p: # labeling initial states
        table_p.update({p: latexname(p,simplify)})
        p1,p2 = p
        if (p1 in Xm1) & (p2 in Xm2):
            Xm_p = Xm_p | frozenset([p])        
    # main stack cycle     
    while S != [] :
        p = S.pop(); # taking a tuple from stack 
        p1, p2 = p   #     
        # evaluating in tha active event set
        for sigma in Gamma_prod(p):
            Q1_n, Q2_n = delta1(p1, sigma), delta2(p2,sigma)
            Q = frozenset((q1,q2) for q1 in Q1_n  for q2 in Q2_n)
            # Q contains the tuples fo reached states
            for q in Q:                 
                if q not in X_p: 
                    # updating states of composite automaton and stack                                                    
                    X_p = X_p | frozenset([q])           
                    S.append(q)
                    q1, q2 = q
                    # updating marked states
                    if (q1 in Xm1) & (q2 in Xm2):
                        Xm_p =Xm_p | frozenset([q]) 
                    #updating latex dictionary
                    table_p.update({q: latexname(q,simplify)})
                transition.append([p,sigma,q])  
            if sigma in self.symDict:
                table_p.update({sigma: self.symDict[sigma]}) 
            elif sigma in other.symDict:
                table_p.update({sigma: other.symDict[sigma]})  


    #Simplify
    if simplify:
        X_p2=()
        for i in X_p:
            X_p2+=(pars(i),)
        X_p = X_p2
        for i in range(len(transition)):
            for j in range(len(transition[i])):
                if type(transition[i][j])==tuple:
                    transition[i][j]=pars(transition[i][j])
        for i in range(len(X0_p)):
            if type(X0_p[i])==tuple:
                X0_p[i]=pars(X0_p[i])
        Xm_p=list(Xm_p)
        for i in range(len(Xm_p)):
            if type(Xm_p[i])==tuple:
                Xm_p[i]=pars(Xm_p[i])
        Xm_p=frozenset(Xm_p)
        keys=list(table_p.keys())
        values=list(table_p.values())
        for i in range(len(keys)):
            if type(keys[i])==tuple:
                keys[i]=pars(keys[i])
        table_p={}
        for i in range(len(values)):
            if values[i].count(")")>0:
                values[i]="".join(values[i].split("("))
                values[i]="".join(values[i].split(")"))
                values[i]="("+values[i]+")"
            table_p[keys[i]]=values[i]
    
    # finally, we build the composite automaton    
    prodnondet = fsa(X_p, Sigma_p, transition, X0_p, Xm_p, table=table_p,
                  Sigobs = Sigobs_p, Sigcon=Sigcon_p,name='untitled')    
    return prodnondet



def productdet(self,other,simplify=True):
    """This function builds an automaton corresponding
    to the product of the input automata.
    Specifically, productdet is optimized for deterministic
    automata and it generates an error if the user tries to
    use it with nondeterministic automata"""    
    from deslab.src.automatadefs import fsa
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
        if x1 in self.symDict:
            x1_name = self.symDict[x1]
        else: 
            x1_name = str(x1)
        if x2 in other.symDict:
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
    # easy names for inner variables
    from deslab.src.comparison import isitempty
    from deslab.src.automatadefs import fsa
 
    # G1 marks an empty language and G2 does not
    if self.empty and not other.empty:
        return fsa()
    # G1 and G2 marks empty languages
    elif self.empty and other.empty:
        return fsa()
    # G1 marks a nonempty language and G2 do.
    elif not self.empty and other.empty:
        return fsa() 
        
    X1, X2 = self.X, other.X
    X01, X02 = self.X0, other.X0
    Xm1, Xm2 = self.Xm, other.Xm 
    Sigma1, Sigma2 = self.Sigma, other.Sigma
    Gamma1, Gamma2 = self.Gamma, other.Gamma 
    delta1, delta2 = self.delta, other.delta   
    Sigobs1,Sigobs2 = self.Sigobs, other.Sigobs         
    Sigcon1,Sigcon2 = self.Sigcon, other.Sigcon  
    # initial values for composite automaton
    p0 = (list(X01)[0], list(X02)[0])
    X0_p = [p0]
    S = [p0]
    X_p = frozenset(X0_p)
    Xm_p = EMPTYSET
    Sigobs_p = Sigobs1|Sigobs2 
    Sigcon_p = Sigcon1|Sigcon2
    Sigma_p = Sigma1 | Sigma2
    transition=[]
    table_p = {}
    table_p.update({p0: latexname(p0,simplify)})
    p01, p02 = p0
    # is initial state marked ?
    if (p01 in Xm1) & (p02 in Xm2):                
        Xm_p = Xm_p | frozenset([p0])     
    # main stack cycle    
    while S != [] :
        p = S.pop(); # taking a tuple
        p1, p2 = p     
        for sigma in Gamma_prod(p):
            q = delta1(p1, sigma), delta2(p2, sigma) # reached state                                 
            if q not in X_p:
                # updating state and stack                                                     
                X_p = X_p | frozenset([q])           
                S.append(q)
                q1, q2 = q
                if (q1 in Xm1) & (q2 in Xm2):
                    # updating marked states
                    Xm_p =Xm_p | frozenset([q])
                # updating labeling table 
                table_p.update({q: latexname(q,simplify)})
            transition.append([p,sigma,q])  
            if sigma in self.symDict:
                table_p.update({sigma: self.symDict[sigma]}) 
            elif sigma in other.symDict:
                table_p.update({sigma: other.symDict[sigma]})  
    # building the composite automaton
    if simplify:
        X_p2=()
        for i in X_p:
            X_p2+=(pars(i),)
        X_p = X_p2
        for i in range(len(transition)):
            for j in range(len(transition[i])):
                if type(transition[i][j])==tuple:
                    transition[i][j]=pars(transition[i][j])
        for i in range(len(X0_p)):
            if type(X0_p[i])==tuple:
                X0_p[i]=pars(X0_p[i])
        Xm_p=list(Xm_p)
        for i in range(len(Xm_p)):
            if type(Xm_p[i])==tuple:
                Xm_p[i]=pars(Xm_p[i])
        Xm_p=frozenset(Xm_p)
        keys=list(table_p.keys())
        values=list(table_p.values())
        for i in range(len(keys)):
            if type(keys[i])==tuple:
                keys[i]=pars(keys[i])
        table_p={}
        for i in range(len(values)):
            if values[i].count(")")>0:
                values[i]="".join(values[i].split("("))
                values[i]="".join(values[i].split(")"))
                values[i]="("+values[i]+")"
            table_p[keys[i]]=values[i]


    proddet = fsa(X_p, Sigma_p, transition, X0_p, Xm_p, table=table_p,
                  Sigobs = Sigobs_p, Sigcon=Sigcon_p,name='untitled')    
    return proddet



def parallel(self, other,simplify=True):  
    """This function builds an automaton corresponding
    to the parallel composition  of the input automata. 
    It selects the adequate function to build the parallel
    composition  automaton depending if the input automata
    are deterministic. If they are, it uses paralleltdet,
    otherwise uses parallelnondet

    When simplify=True, the names of the states are
    simplified to a tuple containing the names of the states
    Example
    -------
    #simplify=False
    state = ((x1,x2),x3)

    #simplify=True
    state = (x1,x2,x3)
    
    """  
    from deslab.src.comparison import isitempty
    from deslab.src.automatadefs import fsa
 
    # G1 marks an empty language and G2 does not
    if self.empty and not other.empty:
        return fsa()
    # G1 and G2 marks empty languages
    elif self.empty and other.empty:
        return fsa()
    # G1 marks a nonempty language and G2 do.
    elif not self.empty and other.empty:
        return fsa() 
    # case deterministic, method is_dfa is used
    if self.is_dfa() & other.is_dfa():
        return paralleldet(self,other,simplify)    
    else:# case non deterministic
        return parallelnondet(self,other,simplify)
    
    
def product(self, other, simplify=True): 
    """This function builds an automaton corresponding
    to the product of the input automata. It selects 
    the adequate function to build the product automaton
    depending if the input automata are deterministic. If
    they are, it uses productdet, otherwise uses
    productnondet

    When simplify=True, the names of the states are
    simplified to a tuple containing the names of the states
    Example
    -------
    #simplify=False
    state = ((x1,x2),x3)

    #simplify=True
    state = (x1,x2,x3)
    
    Example
    -------
    syms('a b c d e f s1 s2 s3 s4')
    table=[(a,'\\alpha'),(b,'\\beta'),(c,'c_c'),(s1,'x_1'),(s2,'x_2'),(s3,'x_3'),(s4,'y_4')]
    X1=[s1,s2,s3,s4] #definindo estados
    E1=[a,b,c] #definindo eventoson
    Transition1=[(s1,a,s2),(s2,b,s3),(s3,c,s4),(s1,b,s4)]
    X0=[s1]
    Xm=[s2,s3,s4]
    H1=fsa(X1,E1,Transition1,X0,Xm,table,name='$G_1$')
    X2=[s1,s2,s3,s4] #definindo estados
    E2=[a,b,c,d,e,f] #definindo eventoson
    Transition2=[(s1,a,s2),(s2,e,s3),(s3,f,s4),(s1,b,s3),(s1,f,s4),(s2,f,s4)]
    X0,Xm=[s1],[s1,s2]
    H2=fsa(X2,E2,Transition2,X0,Xm,table,name='$G_2$')
    P = product(H1,H2)
    draw(P2)
 
    """
    from deslab.src.comparison import isitempty
    from deslab.src.automatadefs import fsa 
    # G1 marks an empty language and G2 does not
    if self.empty and not other.empty:
        return fsa()
    # G1 and G2 marks empty languages
    elif self.empty and other.empty:
        return fsa()
    # G1 marks a nonempty language and G2 do.
    elif not self.empty and other.empty:
        return fsa()     
    # case deterministic, method is_dfa is used      
    if self.is_dfa() & other.is_dfa():
        return productdet(self,other,simplify)
    else:# case non deterministic
        return productnondet(self,other,simplify)

def epsilonobserver(self):
    """ This function calculates a deterministic automaton
    which marks and generates the projection of the marked and
    generated languages of the input automaton G over the alphabet
    Sigma - set([epsilon]). It returns the determinized automaton
    with the same languages as the input automaton"""
    from deslab.src.automatadefs import fsa
    # drawing epsilon event from Sigma
    events_noneps = self.Sigma - EMPTYSTRINGSET
    # calculating observer
    epsobserver = observer(self,events_noneps)
    return epsobserver
    


    
