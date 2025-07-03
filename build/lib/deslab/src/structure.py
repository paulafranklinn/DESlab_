"""
Automaton structure manipulation functions for the automata objects.
The purpose of these methods is to provide low lovel access to an
automaton object.
"""
#    Copyright (C) 2011-2012 by
#    Leonardo Bermeo <lbermeoc@unal.edu.co>
#    Joao Carlos Basilio <basilio@poli.ufrj.br>
#    GNU license.
import networkx as nx
import pandas as pd
from deslab.src.def_const import *
import copy
from deslab.src.exceptions import*
from deslab.src.automatadefs import create_FSA_transdicts
set = frozenset


def addtransition(self, trans):   
    """
    This method provides low level access to the transition
    structure of the input automaton. By using this method, user
    can add a new transition (x,e,y) to the automaton G. 
    The states x, y, and the event e are automatically added 
    if they are not already present in the automaton G.
    The output is directed to a different object than G. 
    
    Example
    --------
    
    syms('a b c d e sf x1 x2 x3 x4 x5 x6 x7')
    S = [a,b,c,d,e,sf]
    X=[x1,x2,x3,x4,x5,x6,x7]
    X0=[x1]
    Xm=[x2,x4]
    T=[(x1,c,x2),(x1,a,x5),(x2,sf,x3),(x3,e,x4),(x4,d,x4),
       (x1,a,x5),(x5,b,x6),(x6,d,x6),(x7,e,x7),(x3,a,x7)]
    tab=[(x1,'x_1'),(x2,'x_2'),(sf,'\\sigma_f')]
    G=fsa(X,S,T,X0,Xm,tab) 
    
    # call it as follows
    G1=G.addtransition([x1,b,x2])
    
    # or alternatively
    G1=addtransition(G,[x1,b,x2])
    """

    # making a local copy of the input automaton               
    auto = self.copy()   
    try: 
        x, e, y = trans
    except:
        raise invalidArgument('transition must be a list of the form [x,e,y] or (x,e,y)')
    
    # Case state y of transition is not already in X
    if y not in self.X:
        auto = auto.addstate(y)
    # Case event e of transition is not already in Sigma    
    if e not in self.Sigma:
        auto.Sigma = auto.Sigma | frozenset([e]) 
        auto.Sigobs = auto.Sigobs | frozenset([e]) 
        auto.Sigcon = auto.Sigcon | frozenset([e]) 
    # Case state x of transition is not already in X                      
    if  x  not in  self.X:
        auto = auto.addstate(x)
        # updating graph structure
        auto.Graph.add_edge(x, y, key=e, label=e) 
        # we create a new entry in deltaDict
        auto.deltaDict[x] = {e:frozenset([y])}
        auto.gammaDict[x]=frozenset([e])
    # if event e is an active event of x               
    elif e in self.Gamma(x): 
        # it could be an existing transition
        if y in self.__delta__(x,e):
            return auto  
        # otherwise it could be a new one
        else:      
            # the output automata becomes NDFA
            auto.Graph.add_edge(x, y, key=e, label=e)
            auto.infoDict['isDFA'] = False  
            # we update de deltaDcit entry        
            X_next = self.deltaDict[x][e]  
            X_next = X_next | frozenset([y])
            auto.deltaDict[x][e] = X_next
    # if e is not an active event of x, but x has active events 
    elif auto.Gamma(x) != EMPTYSET:
        # the new transition doesn't change DFA
        auto.Graph.add_edge(x, y, key=e, label=e)
        # we add the new state to transition function
        auto.deltaDict[x][e] = frozenset([y])    
        # ,... and the new event in active events
        auto.gammaDict[x] = self.gammaDict[x] | frozenset([e])
    # if state x does not have active events          
    elif auto.Gamma(x) == EMPTYSET :
        # we add new entries to DFA
        auto.Graph.add_edge(x, y, key=e, label=e)
        # new transition function entry
        auto.deltaDict[x] = {}
        auto.deltaDict[x].update({e:frozenset([y])})
        # new gamma entry and state set update
        auto.gammaDict[x]=frozenset([e])
        auto.X = auto.X | frozenset([y]) 
    else:
        raise deslabError('unexpected  transition')  
    return auto

def deletetransition(self,trans):  
    """
    This method provides low level access to the transition
    structure of input automaton. By using this method, user
    can delete an existing transition in the automaton G. If
    the transition does not exist in G, the method raises an
    exception.
    
    Example
    -------
    syms('a b c d e sf x1 x2 x3 x4 x5 x6 x7')
    S = [a,b,c,d,e,sf]
    X = [x1,x2,x3,x4,x5,x6,x7]
    X0 = [x1]
    Xm = [x2,x4]
    T= [(x1,c,x2),(x1,a,x5),(x2,sf,x3),(x3,e,x4),(x4,d,x4),
        (x1,a,x5),(x5,b,x6),(x6,d,x6),(x7,e,x7),(x3,a,x7)]
    tab = [(x1,'x_1'),(x2,'x_2'),(sf,'\\sigma_f')]
    G = fsa(X,S,T,X0,Xm,tab) 
    
    # call it as follows:
    G1 = G.deletetransition([x1,b,x2])
    G1 = G.deletetransition((x1,b,x2))
    
    # or alternatively:
    G1 = addtransition(G,[x1,b,x2])
    """
    # local copy of input automaton
   
    auto=self.copy()
    try: 
    # only is valid a triple in list or tuple form
        x, e, y = trans
    except:
        raise invalidArgument('transition must be a list of the form [x,e,y] or (x,e,y)')
    try: 
        # we try to set the edge in graph
        auto.Graph.remove_edge(x, y, key=e)           
    except:
        # case it does not exist
        raise invalidTransition('transition %s is does not exist in G'%(str(trans)))
    
    # testing if e labels two differens transitions at x
    # next set of states reached from x, excluding y
    X_next = auto.deltaDict[x][e] - frozenset([y]) 
    # in this case is e labels a single transition        
    # we update delta   
    if X_next == EMPTYSET:
        del  auto.deltaDict[x][e]
        # we delete an empty entry
        if auto.deltaDict[x] == {}:
            del auto.deltaDict[x]
        # actually e leaves gamma
        active_events = auto.gammaDict[x] - frozenset([e])    
        if active_events == EMPTYSET:
            # no more active events
            del  auto.gammaDict[x] 
        else:
            # there are other active events
            auto.gammaDict[x] = active_events   
    else: 
    # in this case there are other transitions
    # labelled by event e       
        auto.deltaDict[x][e] = X_next    
          
    return auto  



def renamevents(self ,mapping): 
    """
    This function renames the events of the input automaton
    according to the input mapping.
    The input mapping is a list or tuples as follows:
    mapping = [(old_e1, new_ e1),(old_e2,new_e2),...]
     
    or it can be used a dictionary as follows
    mapping = {old_e1: new_e1, old_e2: new_e2,...}. 

    
    Example
    --------
    syms('x1 x12 x123 x0 xnot1 x0new ')
    e1,e2,e3 ='1','2','3'
    X1=[x0,x1,x12,x123,xnot1] 
    S=[e1,e2,e3] 
    T = [(x0,e1,x1),(x0,e2,xnot1),(x0,e3,xnot1),(x1,e1,x1),(x1,e2,x12),(x1,e3,xnot1),
         (x12,e1,x1),(x12,e2,xnot1),(x12,e3,x123),(x123,e1,x1),(x123,e2,xnot1),
         (x123,e3,xnot1),(xnot1,e1,x1),(xnot1,e2,xnot1),(xnot1,e3,xnot1)]    
    X0=[x0]
    Xm=[x123]
    G1=fsa(X1,S,T,X0,Xm,name='$G_1$')
    map=[(e1,a),(e2,b)]
    
    # call it as follows
    G2=G1.renamevents(map)
    
    # or alternatively
    G2=renamevents(G1, map)
    """  
    
    from deslab.src.automatadefs import create_FSA_transdicts
    # we make a copy of input automaton
    auto = self.copy()
    # only dictionaries or lists are allowed
    if isinstance(mapping,list):
        mapping=dict(mapping)
    if not isinstance(mapping,dict):
        raise invalidArgument('Mapping must be a list of relations, string or dictionary')
    # we look for transitions labelled by e    
    for x in auto:
        for e in mapping:
            #updating the symbolic dictionary
            if e in auto.symDict:
                del auto.symDict[e]   
            # finding the next state in each transition         
            if e in auto.Gamma(x):
                # next state set reached from state x
                X_next = auto.deltaDict[x][e]
                for x_n in X_next:   
                    # updating Graph                 
                    auto.Graph.add_edge(x, x_n, key=mapping[e], label=mapping[e])
                    auto.Graph.remove_edge(x,x_n,e)  
                    # updating Gamma dictionary  
                    act_events = auto.gammaDict[x]
                    act_events = act_events | frozenset([mapping[e]])
                    act_events = act_events - frozenset([e])                    
                    auto.gammaDict[x] = act_events
                    # updating delta dictionary
                    auto.deltaDict[x][mapping[e]] = auto.deltaDict[x].pop(e)
    # updating Sigma                
    events_not_ren = auto.Sigma - set(mapping.keys())   
    auto.Sigma = events_not_ren | set(mapping.values())    
        
    # updating observable events
    new_obs = set([mapping[e] for e in auto.Sigobs if e in mapping])
    old_obs = set([e for e in auto.Sigobs if e not in mapping])
    auto.Sigobs =  new_obs | old_obs 

    #updating controllable events
    new_con = set([mapping[e] for e in auto.Sigcon if e in mapping])
    old_con = set([e for e in auto.Sigcon if e not in mapping])
    auto.Sigcon =  new_con | old_con
    return auto


def renamestates(self, mapping):
    """
    This function renames the events of the input automaton
    according to the input mapping.
    The input mapping is a list or tuples as follows:
    mapping = [(old_x1, new_ x1),(old_x2,new_x2),...]
     
    or it can be used a dictionary as follows
    mapping = {old_x1: new_x1, old_x2: new_x2,...}. 
    
    Example
    -------
    
    syms('x1 x12 x123 x0 xnot1 x0new  y1 y3 y12 y123 y0')
    e1,e2,e3 = 1,2,3
    X1=[x0,x1,x12,x123,xnot1] 
    S=[e1,e2,e3] 
    T = [(x0,e1,x1),(x0,e2,xnot1),(x0,e3,xnot1),(x1,e1,x1),(x1,e2,x12),(x1,e3,xnot1),
         (x12,e1,x1),(x12,e2,xnot1),(x12,e3,x123),(x123,e1,x1),(x123,e2,xnot1),
         (x123,e3,xnot1),(xnot1,e1,x1),(xnot1,e2,xnot1),(xnot1,e3,xnot1)]
    X0=[x0]
    Xm=[x123]
    G1=fsa(X1,S,T,X0,Xm,name='$G_1$')
    map=[(x1,y1),(x12,y12),(x123,y123),(xnot1,x0new),(x0,y0)]

    # call it as follows
    G2 = G1.renamestates(map)
    
    # or alternatively
    G2 = renamestates(G1, map)
    
    

    
    
    """ 
 
    # local copy of automaton G
    auto = self.copy()      
    # we revise correctness of input parameters   
    if isinstance(mapping,list):
        mapping=dict(mapping)    
    elif mapping=='number':
        mapping = lexgraph_numbermap(auto) 
    elif mapping=='lex':
        mapping = lexgraph_alphamap(auto)                               
    if not isinstance(mapping,dict):
        raise invalidArgument('Mapping must be a list of relations, string or dictionary')
    
    # updating the names in the  graph structure       
    auto.Graph = nx.relabel_nodes(auto.Graph,mapping)
    # changing elements of the sets
    for x in mapping:
        # updating the symbolic dictionary
        if x in self.symDict:
            del auto.symDict[x] 
        # updating set of states X 
        auto.X = auto.X - frozenset([x]);
        auto.X =auto.X | frozenset([mapping[x]]) 
        # updating set of marked states Xm
        if x in self.Xm:
            auto.Xm = auto.Xm - frozenset([x]);
            auto.Xm =auto.Xm | frozenset([mapping[x]])             
        # updating set of initial states
        if x in self.X0:
            auto.X0 = auto.X0 - frozenset([x]);
            auto.X0 = auto.X0 | frozenset([mapping[x]])
    auto.gammaDict, auto.deltaDict, auto.infoDict = create_FSA_transdicts(auto.Graph,auto.X,auto.Sigma,auto.X0)
    return auto

    


def addevent(self, sigma):
    """
    This method extends the input alphabet of the automaton G1 without
    modification of the transition structure of G1. It is used for
    standarsizing the alphabets in common operations of two input 
    automata. The new events added are defined both controllable and
    observable. Use the method setpar to change the definition of 
    these sets.
    
    
    Example:
    -------
    
    syms('q1 q2 q3 q4 a1 b1 e f')
    table = [(a1,'a_1'),(b1,'b_1'),(q1,'q_1'),(q2,'q_2'),(q3,'q_3'),(q4,'t_1')]
    X = [q1,q2,q3,q4]
    Sigma = [a1,b1]
    X0 = [q1]
    Xm = [q1,q3]
    T =[(q1,a1,q2),(q1,b1,q3),(q2,b1,q3),(q3,b1,q4)]
    G1 = fsa(X,Sigma,T,X0,Xm,table,name='$G_1$ -- example 1')
    Gext= G1.addevents(e)
    or
    Gext= G1.addevents([e,f])
    print Gext.Sigma
    """
    #we make a inner copy of the input automaton
    auto = self.copy()
    if type(sigma)==str:
        newset = frozenset([sigma])
    elif type(sigma) in [list, tuple, frozenset]:
        newset = frozenset(sigma)
    
    # we check for validity of the data type of set passed
    try:         
        auto.Sigma = auto.Sigma | newset
    except:
        raise invalidArgument('event sigma is not valid')
    # ... the new evens are considered observable by default
    auto.Sigobs = auto.Sigobs | newset 
    # ... and also controllable
    auto.Sigcon = auto.Sigcon | newset 
    return auto


def deletevent(self, sigma):
    """  
    This method provides low level access to the transition
    structure of input automaton. By using this method, user
    can delete the sigma event in the alphabet of automaton G.
    All transitions labelled by sigma will be deleted
    
    Example
    -------
    syms('a b c d e sf x1 x2 x3 x4 x5 x6 x7')
    S = [a,b,c,d,e,sf]
    X = [x1,x2,x3,x4,x5,x6,x7]
    X0, Xm  = [x1], [x1]
    T = [(x1,c,x2),(x1,a,x5),(x2,sf,x3),(x3,e,x4),(x4,d,x4),(x1,a,x5),
         (x5,b,x6),(x6,d,x6),(x7,e,x7),(x3,a,x7),(x1,a,x2)]
    G = fsa(X,S,T,X0,Xm) 
    G1=G.deletevent(a)
    draw(G,G1)     
    """
    # it raises an exception if is not possible to delete sigma
    if sigma not in self.Sigma:
        raise invalidArgument('event  %s is not in the input alphabet of G'%(str(sigma)))
        
    # we make a local copy of input automaton
    auto=self.copy()     
    # drawing sigma of the input alphabet
    auto.Sigma = auto.Sigma- frozenset([sigma])
    auto.Sigobs = auto.Sigobs- frozenset([sigma])
    auto.Sigcon = auto.Sigcon- frozenset([sigma]) 
    # we find the transitions labelled by e
    for edge in self.Graph.edges(keys=True):
        x, y, e = edge
        # if any founded, delete it
        if (e == sigma) and (x in auto.gammaDict):
            # removing sigma from graph
            auto.Graph.remove_edge(x, y, key=e)
            act_events = auto.Gamma(x) - frozenset([e])           
            if act_events == EMPTYSET:
                # all transitions of x exhausted
                del auto.gammaDict[x]
                del auto.deltaDict[x]
            else: 
                # still surviving some transitions
                auto.gammaDict[x] = act_events  
                X_next = auto.deltaDict[x][e] - frozenset([y])
                if X_next == EMPTYSET:
                    # if states accessible from x through sigma
                    # are exhausted
                    del auto.deltaDict[x][e]
                else:
                    # still surviving transitions
                    auto.deltaDict[x][e] = X_next                                            
    return auto 















def addstate(self, x, marked = False): 
    """
    This function adds the state x to the set of states of
    automaton G. If the parameter marked is set to True, the new 
    state added is marked, otherwise it will be a nonmarking
    state. By default the parameter marked is False.
    
    Example 
    
    from  deslab import *

    syms('q1 q2 q3 q4 x1 y1 a1 b1 e f')
    X = [q1,q2,q3]
    Sigma = [a1,b1]
    X0 = [q1]
    Xm = [q1,q3]
    T =[(q1,a1,q2),(q1,b1,q3),(q2,b1,q3)]
    G1 = fsa(X,Sigma,T,X0,Xm)
    G2 = G1.addstate(q4)
    draw(G1,G2)
        
    """
    # we make a local copy of input automaton
    auto = self.copy()
    labels = [int(i[1]['label'][1::]) for i in auto.Graph.nodes(data=True)]
    i = max(labels)+1
    try:    
        newstate = frozenset([x])  
        auto.Graph.add_node(x, label='s'+str(i))          
        auto.X = auto.X | newstate
        
    except: 
        raise deslabError('state  %s is not a  valid expression'%(str(x)))
    if marked:
        auto.Xm = auto.Xm | newstate
    return auto

def deletestate(self, x):     
        
    from deslab.src.automatadefs import fsa
    # we make a local copy of input automaton
    auto=self.copy()      
    if x not in self.X:
        # we raises an error for weird state
        raise deslabError('state  %s does not belong to the set of states in input automaton'%(str(x)))
    if len(self) == 1:
        # a one state automata becomes empty
        return fsa()
    # deleting outcoming edges
    if x in auto.gammaDict:
        del auto.gammaDict[x]
        del auto.deltaDict[x]      
    # updating sets
    auto.Graph.remove_node(x)            
    auto.X = auto.X - frozenset([x])
    auto.Xm = auto.Xm - frozenset([x])
    auto.X0 = auto.X0 - frozenset([x])
    # updating incoming edges    
    for edge in self.Graph.in_edges(x,keys=True):
        e, x_p = edge[2], edge[0]   
        # predecessor state has active states
        if x_p  in auto.gammaDict:           
            # event e is one of the active estates          
            if e in auto.deltaDict[x_p]:          
                # next state set for nondeterministic automata,
                # excluding  state x  deleted  
                X_next = auto.deltaDict[x_p][e] - frozenset([x])
                # active events of x_p            
             
                if X_next == EMPTYSET:
                # In the case of no other transitions
                # we erase the transition
                    del auto.deltaDict[x_p][e]
                    # if it is the last transition 
                    if auto.deltaDict[x_p] == {}:
                        del auto.deltaDict[x_p]
                    # in this case there are not other transitions labelled
                    # by e, then update the gamma dictionary, drawing that event
                    active_events = auto.gammaDict[x_p] - frozenset([e])
                    if active_events == EMPTYSET:
                        del  auto.gammaDict[x_p] 
                    else:
                        auto.gammaDict[x_p] = active_events                                        
                else:
                    # in the case of more transitions
                    # labelled by event e, we update delta
                    auto.deltaDict[x_p][e] = X_next             
    
    return auto     

def renametransition(G, rentrans):
    """This method changes the event that labels a specific
    transition. Suppose that automaton G possess a transition
    (x, e_old, y), by using this method the user can set a new 
    event labeling this transition 
    
    Example
    -------
    
    syms('q1 q2 q3 q4 x1 y1 a1 b1 e f')
    X = [q1,q2,q3,q4]
    Sigma = [a1,b1]
    X0 = [q1]
    Xm = [q1,q3]
    T =[(q1,a1,q2),(q1,b1,q3),(q2,b1,q3),(q3,b1,q4)]
    G1 = fsa(X,Sigma,T,X0,Xm)
    Gr=G1.renametransition([q1,(a1,f),q2])
    draw(G1,Gr)    
    """    
    self = G
    try:
        x, pair, y = rentrans 
        oldname, newname =  pair 
    except:
        raise invalidArgument('trans must be a list of the form [x, (old_e, new_e), y]')
    auto = self.deletetransition((x, oldname,y))
    auto = auto.addtransition((x, newname ,y ))
    return auto

def addselfloop(G, x, e):
    
    """
    This method adds a selfloop to the state x which
    is labelled by the event e.
    Example
    -------
    
    syms('q1 q2 q3 q4 x1 y1 a1 b1 e f')
    X = [q1,q2,q3,q4]
    Sigma = [a1,b1]
    X0 = [q1]
    Xm = [q1,q3]
    T =[(q1,a1,q2),(q1,b1,q3),(q2,b1,q3),(q3,b1,q4)]
    G1 = fsa(X,Sigma,T,X0,Xm)
    Gsl=G1.addselfloop(q1,a1)
    draw(G1,Gsl)    
    """  
    auto = G.addtransition([x,e,x])
    return auto
  
def transitions_iter(G):
    """
    This method iterates over the transitions of G.
    
    Example
    -------
    
    from deslab import *
    syms('q1 q2 q3 q4 x1 y1 a1 b1 e f')
    X = [q1,q2,q3,q4]
    Sigma = [a1,b1]
    X0 = [q1]
    Xm = [q1,q3]
    T =[(q1,b1,q1),(q1,a1,q2),(q1,b1,q3),(q2,b1,q3),(q3,b1,q4)]
    G1 = fsa(X,Sigma,T,X0,Xm)
    for i in transitions_iter(G1):
        print i
    """  
    for T in G.Graph.edges(keys=True):
        yield [T[0], T[2], T[1]]
    
  
def transitions(G):        
    """
    This method returns the list of transitions that
    form the structure of automaton G
    
    Example
    -------
    
    from deslab import *
    syms('q1 q2 q3 q4 x1 y1 a1 b1 e f')
    X = [q1,q2,q3,q4]
    Sigma = [a1,b1]
    X0 = [q1]
    Xm = [q1,q3]
    T =[(q1,b1,q1),(q1,a1,q2),(q1,b1,q3),(q2,b1,q3),(q3,b1,q4)]
    G1 = fsa(X,Sigma,T,X0,Xm)
    T = transitions(G1)
    print T 
    """      
    T = [(t[0], t[2], t[1]) for t in G.Graph.edges(keys=True)]
    return T
      
def lexgraph_dfs(G):    
    """Return a list of states whose order is given by the 
    lexicografical depth first search from initial state """
    x0 = list(G.X0)[0]
    events_x0 = sorted(list(G.Gamma(x0)))
    stack = []
    visited = [x0]   
    node = [epsilon]             
    stack = [(x0, iter(events_x0))]              
    while stack:
        x, event_iter = stack[-1]   
        try:
            e = next(event_iter)
            x_n = G.delta(x, e)        
            if x_n  not in visited:                       
                visited.append(x_n)
                events_new = sorted(list(G.Gamma(x_n)))                
                stack.append((x_n,iter(events_new)))             
        except StopIteration:
            stack.pop()                 
    return visited


def lexgraph_alphamap(G):  
    """Return a mapping of states obtained when is performed 
    a lexicographical depth first search in automaton G. The output
    is a dictionary that associates every state x of G1 with the string
    formed with events of the minimal path that conduces to x. 
     """  
    x0 = list(G.X0)[0]
    events_x0 = sorted(list(G.Gamma(x0)))
    stack = []
    visited = [x0]   
    node = [] 
    mapping = {x0: epsilon}           
    stack = [(x0, iter(events_x0))]              
    while stack:
        x, event_iter = stack[-1]   
        try:
            e = next(event_iter)
            x_n = G.delta(x, e)        
            if x_n  not in visited:                       
                node += [e]
                visited.append(x_n)
                mapping.update({x_n: "".join(node)})                
                events_new = sorted(list(G.Gamma(x_n)))                
                stack.append((x_n,iter(events_new))) 
        except StopIteration:
            node = node[0:-1]
            stack.pop()           
    return mapping


def lexgraph_numbermap(G):    
    """Return a mapping of states obtained when is performed 
    a lexicographical depth first search in automaton G. The output
    is a dictionary that associates every state x of G1 with a number 
    between 1 and N (number of states of G1) that represent the order
    of the state x according to the lexicographical depth first search. 
    """      
    x0 = list(G.X0)[0]
    events_x0 = sorted(list(G.Gamma(x0)))
    visited = [x0]   
    node = 0 
    mapping = {x0: node}           
    stack = [(x0, iter(events_x0))]   
    while stack:
        x, event_iter = stack[-1]   
        try:
            e = next(event_iter)
            x_n = G.delta(x, e)
            if x_n  not in visited:                         
                node = node + 1
                visited.append(x_n)
                mapping.update({x_n: node})                
                events_new = sorted(list(G.Gamma(x_n)))                
                stack.append((x_n,iter(events_new))) 
        except StopIteration:
            stack.pop()  
    return mapping

def size(G):
    """return the number of transitions of automaton G"""
    return G.Graph.size()    
    

def mtable(G,states = [], events = [], gen_csv = False, csv_name = 'new'):
    """
    A DataFrame representing the transition table with states as rows and
    events as columns. The table cells contain the resulting state after 
    applying the event to the current state, or None if no transition exists.

            | Event1  | Event2

    State1  | State2  | None

    State2  | None    | State1

    If you want a specific order in the table, provide it as an argument
    
    Also, if you want to generate a csv from your table, input gen_csv = True,
    you can choose the name by given a string to csv_name variable.
    
    -------
    Example
    
    
    X = [q1,q2]
    Sigma = [a1,b1]
    X0 = [q1]
    Xm = [q1,q2]
    T =[(q1,b1,q1),(q1,a1,q2),(q2,b1,q1)]
    G1 = fsa(X,Sigma,T,X0,Xm)
    table = mtable(G1)
    print(table)
    
            | a1      | b1

    q1      | q2      | q1

    q2      | None    | q1    
    """
    if not states:
        states = list(G.X)
    if not events:
        events = list(G.Sigma)
    
    trans = list(G.transitions())

    table = pd.DataFrame(index=states, columns=events) 
    
    for begin, event, end in trans:
        table.loc[begin, event] = end
    
    table = table.where(pd.notnull(table), None)
   
    if gen_csv:
        name = csv_name + '.csv'
        table.to_csv(name)
   
    return table
