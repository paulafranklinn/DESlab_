"""Base class for automata."""
#    Copyright (C) 2011-2012 by
#    Leonardo Bermeo <lbermeoc@unal.edu.co>
#    Joao Carlos Basilio <basilio@poli.ufrj.br>
#    BSD license.
from deslab.src.exceptions import *
import  networkx as nx
import copy
from deslab.src.def_const import *
from deslab.graphics.drawing import graphic

__author__ = """\n""".join(['Leonardo Bermeo (lbermeoc@unal.edu.co)',
                            'Joao Carlos Basilio (basilio@poli.ufrj.br)'])




class fsa: 
    """
    This instruction defines an automaton object inside DESlab.
        
    Parameters
    ----------
    
    X :  list containing the symbols corresponding to states
    Sigma : list of symbols of the input alphabet  
    T : list of transitions ordered as tuples. 
    X0 : list of initial states
    Xm : list of marked states
    table : list of association tuples between symbols and
            latex labels.
    Sigobs : list of observable states
    Sigcon : list of controllable states
    name : a string that represents the name of the system
    type : a string that represents the type of the system
    
    See Also
    --------
    Monoid       
  
    """
                           
    def __init__(self, X=[], Sigma=[], transition=[], X0 =[], Xm=[], 
                 table=[], Sigobs = ALL_EVENTS, Sigcon = ALL_EVENTS, name='untitled', graphic=graphic()):
        
        """
        this is the constructor of the automaton class""" 
        if (X == []) | (X0 == []) | (Sigma == []) : # in the case of empty automaton 
            self.X0 = EMPTYSET
            self.empty=True
            self.Sigma = EMPTYSET            
            self.Xm = EMPTYSET
            self.X = EMPTYSET
            self.name='Empty Automaton'
            self.graphic = graphic
            return
        if (X == EMPTYSET) | (X0 == EMPTYSET) | (Sigma == EMPTYSET) : # in the case of empty automaton 
            self.X0 = EMPTYSET
            self.empty=True
            self.Sigma = EMPTYSET            
            self.Xm = EMPTYSET
            self.X = EMPTYSET
            self.name='Empty Automaton'
            self.graphic = graphic
            return
        # in the case of nonempty automaton 
        self.X = frozenset(X) # we use frozensets instead of sets 
        self.Sigma = frozenset(Sigma)               
        self.Xm = frozenset(Xm)
        self.X0 = verify_fsa_definition(self.X,self.Sigma,self.Xm,X0) # verification of parameters passed      
        self.name = name  # external properties
        self.empty = False 
        
        # controlable and uncontrollable events
        
        if Sigobs == ALL_EVENTS: 
            self.Sigobs = frozenset(Sigma) # default controllable events is Sigma
        else :
            self.Sigobs = frozenset(Sigobs) & self.Sigma    
        if epsilon in Sigma:
            self.Sigobs = self.Sigobs - EMPTYSTRINGSET           
        if Sigcon == ALL_EVENTS: 
            self.Sigcon = frozenset(Sigma) # default controllable events is Sigma
        else :
            self.Sigcon = frozenset(Sigcon) & self.Sigma            
        # constructors of Graph class and 
        self.Graph  = create_graph(transition,self.X,self.Sigma)   
        # delta and Gamma data structures        
        self.gammaDict,self.deltaDict,self.infoDict = create_FSA_transdicts(self.Graph,self.X,self.Sigma,self.X0)
        # symbolic latex labels
        self.symDict = create_table(self.X,self.Sigma, table)   # symbol dictionary        
        # Basic graphic properties, 
        # they should be replaced by an entire class in next version         
        self.graphic = graphic       
        return
       
    def __iter__(self):
        """Iterate over the states. Use the expression 'for n in G'.

        Returns
        -------
        niter : iterator
            An iterator over all states in the automaton.

        Examples
        --------
        """
        return iter(self.X)
    
    def __len__(self):
        """Return the number of states of automaton G.
           Use the expression 'len(G)'.
        """
        return len(self.X)
    
    def info(self):        
        isDFA,hasEpsilon,nonDetTrans = self.infoDict['isDFA'],self.infoDict['hasEpsilon'], self.infoDict['nonDetTrans']
        return isDFA,hasEpsilon,nonDetTrans
        
    def tmx(self,table=None):
        """
        Return the transition matrix of the automaton.
        If table = 'table', prints the table instead

        Example
        -------
        
        syms('a b c d e sf x1 x2 x3 x4 x5 x6 x7')
        sigma=[a,b,c,d,e,sf]
        sigmaf=set([sf])
        X=[x1,x2,x3,x4,x5,x6,x7]
        X0=[x1]
        Xm=[x2,x4]
        T=[(x1,c,x2),(x1,a,x5),(x2,sf,x3),(x3,e,x4),
        (x4,d,x4),(x1,a,x5),(x5,b,x6),(x6,d,x6),(x7,e,x7),(x3,a,x7)]
        table=[(x1,'x_1'),(x2,'x_2'),(sf,'\\sigma_f')]
        G=fsa(X,sigma,T,X0,Xm,table) 

        G.tmx()
        or
        G.tmx('table')
        """
        
        if table=="table":
            text = ''
            m = self.tmx()
            text+="\nTransition Matrix:\n"
            for i in m:
                    text+="\n"
                    for j in i:
                            text+= str(j)+'\t'
            text+="\n"
            print(text)
            return
        
        x = [""]+list(self.X)
        sigma = [""]+list(self.Sigma | self.Sigobs)
        trans = self.transitions()
        m = [['N/D' for j in range(len(sigma))] for i in range(len(x))]
        m[0] = sigma

        for i in range(len(x)):
            m[i][0]=x[i]

        for tup in trans:
            i = x.index(tup[0])
            j = sigma.index(tup[1])
            if m[i][j]!='N/D':
                m[i][j] = m[i][j]+', '+str(tup[2])
            else:
                m[i][j] = str(tup[2])
        return m
    
    def deletetransition(self,trans): 
        """
        This method provides low level access to the transition
        structure of input automaton. By using this method, user
        can delete an existent transition in the automaton G. If
        the transition does not exist in G, the method raises an
        exception.
        
        Example
        -------
        
        syms('a b c d e sf x1 x2 x3 x4 x5 x6 x7')
        sigma=[a,b,c,d,e,sf]
        sigmaf=set([sf])
        X=[x1,x2,x3,x4,x5,x6,x7]
        X0=[x1]
        Xm=[x2,x4]
        T=[(x1,c,x2),(x1,a,x5),(x2,sf,x3),(x3,e,x4),
                    (x4,d,x4),(x1,a,x5),(x5,b,x6),(x6,d,x6),(x7,e,x7),(x3,a,x7)]
        table=[(x1,'x_1'),(x2,'x_2'),(sf,'\\sigma_f')]
        G=fsa(X,sigma,T,X0,Xm,table) 
        G2=G.deletetransition([x7,e,x7])
        draw(G,G2)    
        """         
        from deslab.src.structure import deletetransition
        auto = deletetransition(self, trans)
        return auto
    
    
    
    def addtransition(self, trans):
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
        from deslab.src.structure import addtransition
        auto = addtransition(self, trans)
        return auto
    
    
    
    
    def addevent(self, sigma):
        """
        This method extends the input alphabet of the automaton G1
        without modification of the transition structure of G1. It is used
        for standarsizing the alphabets in common operations of two input 
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
        G11= G1.addevents(e)
        print G11.Sigma
        """
        from deslab.src.structure import addevent
        auto = addevent(self, sigma)
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
        
        from deslab.src.structure import deletevent
        return deletevent(self, sigma)
       
    def addstate(self, x, marked=False):
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
        from deslab.src.structure import addstate        
        return addstate(self, x, marked)
    
    def deletestate(self, x):
        from deslab.src.structure import deletestate        
        return deletestate(self, x)
    

        
    def renametransition(self, trans):
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
        from deslab.src.structure import renametransition        
        return renametransition(self, trans) 
    
    
    def addselfloop(self, x, e):          
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
        
        from deslab.src.structure import addselfloop        
        return addselfloop(self, x, e)     


    def transitions(self):        
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
        T = G1.transitions(G1)
        print T 
        """      
        from deslab.src.structure import transitions        
        return transitions(self)
    
    def transitions_iter(self):        
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
        for i in transitions_iter(G1):
            print i
 
        """      
        from deslab.src.structure import transitions_iter        
        return transitions_iter(self)


    def setpar(self, **param):
        auto = self.copy()
        for key in  param.keys():
            if key =='Xm':
                auto.Xm = frozenset(param[key])
            elif key == 'Sigma':
                sigdiff = auto.Sigma - set(param[key])
                if sigdiff:
                    for event in list(sigdiff):
                        auto = auto.deletevent(event)
                else:
                    sigdiff = set(param[key]) - auto.Sigma
                    auto = auto.addevent(sigdiff)
            elif key =='X0':
                auto.X0 = verify_fsa_definition(auto.X,auto.Sigma,auto.Xm,param[key])
                if len(auto.X0)>1:
                    auto.infoDict['isDFA'] = False
                    auto.infoDict['nonDetTrans'] = True
            elif key =='name':
                auto.name = param[key]
            elif key == 'Sigcon':
                auto.Sigcon = frozenset(param[key])
            elif key == 'Sigobs':
                auto.Sigobs = frozenset(param[key])
            elif key == 'table':
                auto.symDict.update(dict(param[key]))
            else:
                raise invalidArgument('key %s does not exist'%key)
        return auto
    
    def setgraphic(self, style = 'normal', program = 'dot', ranksep = 0.25, nodesep= 0.25, direction = 'LR',
                 FillColor=('plantfill',76), LineColor= ('plantline',85)):
        self.graphic = graphic(style, program, ranksep, nodesep, direction, FillColor, LineColor)
        return
    
            
    def __str__(self):
        text = "finite state machine: %s \n \
        number of states: %s \n \
        number of events: %s \n \
        number of transitions %s \n \
        "%(self.name,len(self.X),len(self.Sigma),self.Graph.size())
        return text

    
    def __delta__(self,state,event):                               
        """This is a special transition function  used in several methods. It keeps
         regularity on data types, and the output state is always a set"""
        if event in self.Sigma: # valid event            
            if event in self.Gamma(state): # factible event
                x_next = self.deltaDict[state][event] 
                return x_next
            else: return EMPTYSET # unfactible event
        else:  # event undefined for the automaton                       
            raise eventMembershipError('Event %s does not belong to Sigma'%(event))
        
    
    def unobsreach(self,Q,Sigma_o = []):               
        """ This function calculates the unobservable reach set of the set Q
        of states, Sigma_o is the set of observable events """                          
        if isinstance(Q,str):
            Q = frozenset([Q]) 
        elif isinstance(Q,list) | isinstance(Q,set):
            Q = frozenset(Q)
        # initial parameters                               
        S = list(Q) # Stack empty
        UR = Q # initial reached set 
        Sigma, Sigobs = self.Sigma, self.Sigobs
        # determinig unobservable events
        if (Sigma_o == []) :
            Sigma_uobs = Sigma - Sigobs
        elif isinstance(Sigma_o,list) | isinstance(Sigma_o,set):
            Sigma_uobs = Sigma - frozenset(Sigma_o)
        elif isinstance(Sigma_o,frozenset):
            Sigma_uobs = Sigma - Sigma_o
        else: 
            raise invalidArgument('Observed set must be of type set, frozenset or list')
        # main cycle of the search algorithm
        while S != []:
            x_reached=S.pop()
            Q_next = EMPTYSET
            act_events = self.Gamma(x_reached) & Sigma_uobs                          
            for sigma in act_events: 
                # it uses the __delta__ function which uses sets for every automaton               
                Q_next = Q_next | self.__delta__(x_reached,sigma)
            for q in Q_next:             
                if q not in UR:
                    UR = UR | frozenset([q])
                    S.append(q) # stack update 
        return UR 
          
    def deltaobs(self,x,e,Sigma_o=[]):      
        """ This function calculates the observable delta function of
        the state of states, Sigma_o is the set of observable events """  
                                         
        Q_reached = self.unobsreach(x,Sigma_o) 
        Q_next = EMPTYSET                       
        for q_r in Q_reached:
            if e in self.Gamma(q_r):
                Q_next = Q_next | self.__delta__(q_r,e)    
        UR = self.unobsreach(Q_next,Sigma_o)          
        return UR   
    
    
    def delta(self,state,event):    
        """ This is de definition of the delta transition function
        for deterministic and non deterministic automata.
        """
        if self.is_dfa():  
        # in the case of DFA we return a symbol or None
            if event in self.Sigma: # the event is valid
                if event in self.Gamma(state): # the event is factible
                    x_next = self.deltaDict[state][event]
                    x_next = list(x_next)[0]
                    return x_next 
                else: return None # the event is unfactible
            else:
                raise eventMembershipError('Event %s does not belong to Sigma'%(event))
        else: 
        # in the case of NDFA we return set or emptyset
            if event in self.Sigma: # the event is valid
                if event in self.Gamma(state):
                    x_next = self.deltaDict[state][event]
                    #x_next = list(x_next)
                    if len(x_next)==1:
                        return list(x_next)[0]
                    return x_next
                else: return EMPTYSET  # the event is unfactible
            else:
                raise eventMembershipError('Event %s does not belong to Sigma'%(event))
    
    def Gamma(self,state):  
        """ This is the active event function"""
        if not (state in self.X):
            raise stateError('State %s is not defined for this automaton'%str(state))
        if len(self.Graph.out_edges(state))==0:
            return frozenset([])               
        active_events = self.gammaDict[state]       
        return active_events
            
    def renamevents(self,mapping):
        from deslab.src.structure import renamevents
        return renamevents(self,mapping)
    
    def renamestates(self, mapping):
        """
        This function allows renaming of states of the input automaton.
        The input mapin is a list of pair of the form
        mapping = [(oldstate1,newstate1),(oldstate2,newstate2)...]    
        """
        from deslab.src.structure import renamestates
        return renamestates(self,mapping)

    """
    ====================================================================
                       Overloaded functions of automata class
    ====================================================================

    """
    
    def copy(self):            
        auto = copy.deepcopy(self)        
        return auto
    
    def __and__(self,other):
        """This function builds an automaton corresponding
        to the product of the input automata. It selects 
        the adequate function to build the product automaton
        depending if the input automata are deterministic. If
        they are, it uses productdet, otherwise uses
        productnondet   
        
        Example
        -------
        
        syms('a b c d e f s1 s2 s3 s4')
        table=[(a,'\\alpha'),(b,'\\beta'),(c,'c_c'),(s1,'x_1'),(s2,'x_2'),(s3,'x_3'),(s4,'y_4')]
        X1=[s1,s2,s3,s4] #definindo estados
        E1=[a,b,c] #definindo eventoson
        T1=[(s1,a,s2),(s2,b,s3),(s3,c,s4),(s1,b,s4)]
        X0=[s1]
        Xm=[s2,s3,s4]
        H1=fsa(X1,E1,T1,X0,Xm,table,name='$G_1$')
        X2=[s1,s2,s3,s4] #definindo estados
        E2=[a,b,c,d,e,f] #definindo eventoson
        T2=[(s1,a,s2),(s2,e,s3),(s3,f,s4),(s1,b,s3),(s1,f,s4),(s2,f,s4)]
        X0,Xm=[s1],[s1,s2]
        H2=fsa(X2,E2,T2,X0,Xm,table,name='$G_2$')
        P = H2&H1
        draw(P)
        """
        from deslab.src.algorithms import product
        prod = product(self,other)
        return prod
    
    def __add__(self, other):
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
        U= G1+G2
        draw(G1,G2,U)        
        """    
        from deslab.src.algorithms import union
        unionaut = union(self,other)
        return unionaut
    
    def __or__(self, other):
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
        U= G1|G2
        draw(G1,G2,U)        
        """    
        from deslab.src.algorithms import union
        unionaut = union(self,other)
        return unionaut
       
    def __floordiv__(self, other):
        """This function builds an automaton corresponding
        to the parallel composition  of the input automata. 
        It selects the adequate function to build the parallel
        composition  automaton depending if the input automata
        are deterministic. If they are, it uses paralleltdet,
        otherwise uses parallelnondet  
        
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
        P = H1//H2
        draw(P)
        """  
        from deslab.src.algorithms import parallel
        par = parallel(self, other)
        return par
    
    def __mul__(self, other):
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
        Xm2 = [q3]
        Sigma2= [a1,b1, f]
        T2 = [(q1,a1,q2),(q2,b1,q3)]
        G2 = fsa(X,Sigma,T2,X0,Xm2,table,name='$G_2$')
        C = G1*G2
        draw(G1,G2,C)
        """     
        
        from deslab.src.algorithms import concatenation
        concat = concatenation(self, other)
        return concat
        
    def __truediv__(self, other):
        """ This function calculates an automaton that marks
        the language Lm(G1)/Lm(G2), where G1 and G2 are the input
        automata. The resulting automaton is a deterministic automaton
        that marks the quotient of the languages marked by the input
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
        Xm2 = [q3]
        Sigma2= [a1,b1, f]
        T2 = [(q1,a1,q2),(q2,b1,q3)]
        G2 = fsa(X,Sigma,T2,X0,Xm2,table,name='$G_2$')
        C = G1/G2
        draw(G1,G2,C)
        """     
        
        from deslab.src.algorithms import langquotient
        quotient = langquotient(self, other)
        return quotient
    
    def __invert__(self):
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
        from deslab.src.algorithms import complement      
        comp = complement(self)
        return comp
    
    def __sub__(self,other):
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
        D= G2-G1
        draw(G1,G2,D)
        """    
        from deslab.src.algorithms import langdiff
        auto = langdiff(self,other)
        return auto
    
    def __le__(self,other):
        """
        This funtion determines whether or not the language marked by
        automaton G1 is contained in the language marked by
        automaton G2, i.e. if Lm(G1) <= Lm(G2).     
        Example
            
        
        Example
        -------
        
        syms('q1 q2 q3 q4 a1 b1 e f')
        table = [(a1,'a_1'),(b1,'b_1'),(q1,'q_1'),(q2,'q_2'),(q3,'q_3')]
        X = [q1,q2,q3,q4]
        Sigma = [a1,b1,e]
        X0 = [q1]
        Xm = [q1,q3]
        T =[(q1,a1,q2),(q1,b1,q3),(q2,b1,q3),(q3,b1,q4)]
        G1 = fsa(X,Sigma,T,X0,Xm,table,name='$G_1$ -- example 1')
        Xm2 = [q1,q2,q3]
        Sigma2= [a1,b1, f]
        T2 = [(q1,a1,q2),(q1,b1,q2),(q2,b1,q3),(q3,b1,q4)]
        G2 = fsa(X,Sigma2,T2,X0,Xm2,table,name='$G_2$')
        print G1<=G2
        draw(G1,G2)
        """  
        from deslab.src.comparison import issublanguage
        test_inclusion = issublanguage(self,other) 
        return test_inclusion 
    
    
    def __ge__(self, other):    
        """ This funtion determines whether or not the language marked by
        automaton G2 is contained in the language marked by
        automaton G1, i.e. if Lm(G1) >= Lm(G2). 
        
        Example
        -------
        
        syms('q1 q2 q3 q4 a1 b1 e f')
        table = [(a1,'a_1'),(b1,'b_1'),(q1,'q_1'),(q2,'q_2'),(q3,'q_3')]
        X = [q1,q2,q3,q4]
        Sigma = [a1,b1,e]
        X0 = [q1]
        Xm = [q1,q3]
        T =[(q1,a1,q2),(q1,b1,q3),(q2,b1,q3),(q3,b1,q4)]
        G1 = fsa(X,Sigma,T,X0,Xm,table,name='$G_1$ -- example 1')
        Xm2 = [q1,q2,q3]
        Sigma2= [a1,b1, f]
        T2 = [(q1,a1,q2),(q1,b1,q2),(q2,b1,q3),(q3,b1,q4)]
        G2 = fsa(X,Sigma2,T2,X0,Xm2,table,name='$G_2$')
        print G1>=G2
        """      
        from deslab.src.comparison import issublanguage
        test_inclusion = issublanguage(other,self) 
        return test_inclusion 
    
  
    def __eq__(self,other):
        from deslab.src.comparison import are_automataequal
        test_equal = are_automataequal(self,other)
        return test_equal
    
    def __ne__(self, other):
        from deslab.src.comparison import are_automataequal
        test_equal = are_automataequal(self,other)
        return not test_equal  
   
    def is_dfa(self):
        return self.infoDict['isDFA']    
    def has_epsilon(self):
        return self.infoDict['hasEpsilon']
    def has_nondetrans(self):
        return self.infoDict['nonDetTrans']
    
    

 
""" The following function are the constructors of many data structures""" 
    
    
def create_graph(T,X,Sigma):
    
    """ This function creates a graph containing, 
    the names provided by user for events, states
    """
    # we use a multidigraph for supporting the
    # transition structure of G
    Graph = nx.MultiDiGraph()        
    i=0
    # this is the set of events in transitions
    Siguns = Sigma # Sigsen = Sigma sensed 
    for state in X:
        Graph.add_node(state,label='s'+str(i))
        i+=1        
    # we add the transitions to the graph   
    for edge in T : 
        x, e, y = edge
         # we revise correctnesss of T passed      
        if x  not in X : # invalid state 
            raise stateMembershipError("State %s in transition %s does not belong to the state set X"%(x, edge))
        if y  not in X : #invalid state
            raise stateMembershipError("State %s in transition %s does not belong to the state set X"%(y, edge))
        if e  not in Sigma :
            raise eventMembershipError("Event %s in transition %s does not belong to the event E"%(e, edge))
        # we add to the graph
        if  not Graph.has_edge(x, y, e) :
            Graph.add_edge(x, y, key=e, label=e)
            # we count the events used      
            
    return Graph



def create_FSA_transdicts(Graph,X,Sigma,X0):
    """ This is the constructor for Gamma function of automaton. 
        It creates a dictionary for the gamma and delta functions of the automaton
    """    
    infoDict={'nonDetTrans': False, 'isDFA': True, 'hasEpsilon': False}
    if len(X0)>1: 
        infoDict['isDFA'] = False
        infoDict['nonDetTrans'] = True 
    
    def determineGamma(state,fsa_features):
        """Creates subdictionaries of the Gamma and delta functions"""
        activeEvents = []
        nodeDict={}       
        #for transition in Graph.out_edges_iter(state,keys=True):
        for transition in Graph.out_edges(state,keys=True):
            event=transition[2]
            nextState=transition[1] 
            if  event==epsilon :    
                fsa_features['isDFA']=False
                fsa_features['hasEpsilon']=True                                                  
            if event in activeEvents :
                temp = nodeDict[event]
                if isinstance(temp,str):
                    temp = frozenset([temp])
                nodeDict[event] = frozenset([nextState]) | temp
                fsa_features['isDFA']=False
                fsa_features['nonDetTrans']=True   
            else :
                nodeDict.update({event:frozenset([nextState])})
            activeEvents.append(event)                                          
        return activeEvents, nodeDict, fsa_features
    gammaDict={}
    deltaDict={}   
    for state in Graph.__iter__():#nodes_iter():
        gammaAct,nodeDict,infoDict=determineGamma(state,infoDict)  
        if gammaAct != []:
            gammaDict[state]= frozenset(gammaAct) 
        if nodeDict != {}:
            deltaDict[state]=nodeDict                 
    return gammaDict,deltaDict,infoDict


def verify_fsa_definition(X,Sigma,Xm,X0):
    # verifying marked states    
    Xm_wrong =  Xm - X    
    if Xm_wrong != EMPTYSET :
        raise markedSetError("States of %s of Xm do not belong to the state set X"%(Xm_wrong))
    # verifying initial states 
    if isinstance(X0,list) | isinstance(X0,set): #type set
        X0 = frozenset(X0)
    elif isinstance(X0,str) | isinstance(X0,tuple): # other types
        X0=frozenset([X0])
    if not (X0 <= X) : # testing for inclusion
        raise initialStateError('Initial state set %s is not contained in state set X'%(X0))        
    return X0

def create_table(X, Sigma, table): 
    """
    this function cleans the current table of automata.
    It deletes all the symbols that does not belong to
    Sigma or X components of the current automata 
    """

    
    if table == []:  # no table was provided
        sym_dict = {}        
    elif isinstance(table,list):  # a list was provided by the user
        sym_dict = dict(table)
        keys = list(sym_dict.keys())
        for key in keys:
            if key not in (X | Sigma):
                del sym_dict[key]                                    
        return sym_dict
    
    elif isinstance(table,dict): # a dictionary was provided by the user or program call
        sym_dict = table.copy()
        keys = list(sym_dict.keys())
        for key in keys:
            if key not in (X | Sigma):
                del sym_dict[key]   
                    
    else:  # no more options
        raise invalidArgument('table must be a list of tuples or a dictionary')
    sym_dict.update({epsilon:r'\varepsilon'})        
    return sym_dict   








