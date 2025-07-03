#    Copyright (C) 2011-2012 by
#    Leonardo Bermeo <lbermeoc@unal.edu.co>
#    Joao Carlos Basilio <basilio@poli.ufrj.br>
#    BSD license.
from deslab import *
from networkx import selfloop_edges

def diagnoser(G,failevent,ret="GD"):
    X = ['N','Y']
    Sigma = [failevent]
    X0 = ['N']
    Xm = []
    T =[('N',failevent,'Y'), ('Y',failevent,'Y')]
    Al = fsa(X,Sigma,T,X0,Xm,name='$A_l$',Sigobs=[])

    if ret=="GL":
        gl=G//Al
        gl.graphic = graphic('observer')
        return gl
    else:
        return observer(G//Al,G.Sigobs)

def simplify(G):
    g=G.copy()
    mapping=[]
    for state in g.X:
        tex="".join(str(i) for i in state)
        tex=tex.replace('(',"")
        tex=tex.replace(')',"")
        tex=tex.replace(',',"")
        tex=tex.replace(' ',"")
        tex=tex.replace("'","")
        mapping.append((state,tex))
    return g.renamestates(mapping)

def Gscc(G,failevent,sigmas=[]):
    
    """Define Gscc as Gd//Gl
    If len(sigmas) > 1, returns Gscci = Gdi//...//Gl, i = len(sigmas)

    Example:
    -------
    syms('0 1 2 3 4 5 6 c a ad p b f cl e d')
    X = [0,1,2,3,4,5,6]
    Sigma = [c,a,ad,p,b,f,cl,e,d]
    X0 = [0]
    Xm = []
    T =[(0,p,1),(1,ad,2),(1,a,3),(2,cl,4),(3,b,5),(4,f,3),(4,c,6),(5,d,0),(6,e,0)]
    G = fsa(X,Sigma,T,X0,Xm,name='$G2$',Sigobs=[c,a,ad,p,b,cl,e,d])

    Gscc0 = Gscc(G2,f)
    Gscc1 = Gscc(G2,f,[a,ad,p,cl])
    Gscc2 = Gscc(G2,f,[[a,ad,p,cl],[a,ad,p,cl,d]])
    """
    
    Gl = diagnoser(G,failevent,"GL")
    Gl = simplify(Gl)
    Gdi=[]
    if sigmas==[]:
        sigmas=[list(G.Sigobs)]
    elif type(sigmas[0])!=list:
        sigmas=[sigmas]
    for sigma in sigmas:
        g=G.copy()
        g=g.setpar(Sigobs=sigma)
        Gd = diagnoser(g,failevent)
        Gdi.append(simplify(Gd))
    Gd=Gdi[0]
    for gdi in Gdi[1:]:
        Gd=Gd//gdi
    Gscc = Gd//Gl
    sigobs = Gscc.Sigobs - set([failevent])
    Gscc = Gscc.setpar(Sigobs=sigobs)
    Gscc.graphic = graphic('observer')
    return Gscc

def is_diagnosable(G,failevent,sigmas=[],method=''):
    """Returns the (co)diagnosability of G for the given observable sigmas in the choosen method.
    If sigmas=[], the function uses G.Sigobs as sigmas.

    Example:
    --------
    syms('0 1 2 3 4 5 6 a b c f u')
    X = [0,1,2,3,4,5,6]
    Sigma = [a,b,c,f,u]
    X0 = [0]
    Xm = []
    T =[(0,a,1),(1,c,2),(1,b,2),(2,a,2),(2,c,2),(1,f,3),(3,b,4),(4,c,5),(5,a,6),(6,u,6)]
    G = fsa(X,Sigma,T,X0,Xm,name='$g$',Sigobs=[a,b,c])
    
    is_diagnosable(G,f,[[a,b],[a,c]],'Gscc')
    is_diagnosable(G,f,[[a,b],[a,c]],'Gv')
    """
    
    if method == 'Gscc':
        #computes the Gscc automaton of G
        G = Gscc(G,failevent,sigmas)
        
        def N_Y(states):
            """Function to detect the occurence of states 'N' and 'Y'
            in the automaton"""
            #List that holds the information of diagnosability of each
            #Sigobs_i in Gscc's states
            checklist = [0]*(len(states)-1)
            if states[-1].count('Y')>0:
                for i in range(len(states[0:-1])):
                    n = states[i].count('N')
                    if n!=0:
                        checklist[i]+=1
            return checklist

        #Search for the strong components of the automaton
        for comp in strconncomps(G):
            # matrix with the size of the number of Gdi's used in Gscc (Gd1//Gd2//...//Gl)
            diag=[0]*(len(list(G.X)[0])-1)
            if len(comp)>1 or list(comp)[0] in selfloop_edges(G.Graph):
                for state in comp:
                    check = N_Y(state)
                    diag = [diag[i]+check[i] for i in range(len(diag))]
                    if diag.count(0)==0:
                        return False 
        return True

    elif method == 'Gv':
        #computes the Gv automaton of G and 
        sigma = G.Sigma
        G= Gv(G,failevent,sigmas)
        #Search for the strong components of the automaton
        for comp in strconncomps(G):
            for state in comp:
                if ( len(comp)>1 or state in selfloop_edges(G.Graph) )and state[-1].count('Y')>0:
                     events = list(G.Gamma(state) & sigma)
                     for event in events:
                         if G.delta(state,event) in comp:
                             return False

                    
        return True



        
    else:
        print("""Choose a method:\n'Gscc'\n or \n 'Gv'""")
                

def Gv(G,failevent,sigmas=[]):
    """Returns Gv for the given observable sigmas.
    If sigmas is empty, the function uses G.Sigobs as sigmas

    Example:
    --------
    syms('0 1 2 3 4 5 6 a b c f u')
    X = [0,1,2,3,4,5,6]
    Sigma = [a,b,c,f,u]
    X0 = [0]
    Xm = []
    T =[(0,a,1),(1,c,2),(1,b,2),(2,a,2),(2,c,2),(1,f,3),(3,b,4),(4,c,5),(5,a,6),(6,u,6)]
    g = fsa(X,Sigma,T,X0,Xm,name='$g$',Sigobs=[a,b,c])

    G_v = Gv(g,f)
    G_v = Gv(g,f,[a,b])
    G_v = Gv(g,f,[[a,b]])
    G_v = Gv(g,f,[[a,b],[a,c]])
    """
    
    def Ri(sigma,sigmaOi,i):
        #Rename events in sigma that aren't in sigmaoi
        table=[]
        mapping=[]
        for event in sigma:
            if event not in sigmaOi:
                if "_" in event:
                    index = event.index('_')
                    mapping.append((event,event+"R%s"%(i)))
                    table.append(("%sR%s"%(event,i),"%s{%s_{R_%s}}"%(event[0:index],event[index+1:],i)))
                else:
                    mapping.append((event,event+"R%s"%(i)))
                    table.append(("%sR%s"%(event,i),"%s_{R_%s}"%(event,i)))
        return table,mapping

    if sigmas==[]:
        sigmas=[list(G.Sigobs)]
    elif type(sigmas[0])!=list:
        sigmas=[sigmas]
        
    SIGMAn = G.Sigma-set([failevent])
    trans=[]
    #creates self-loops with all events of SIGMAn
    for event in SIGMAn:
        trans.append(('N',event,'N'))
        
    #define the automaton An, a single event automaton
    An = fsa(['N'],SIGMAn,trans,['N'],[],name="$A_n$")
    Gn=G&An
    Gn = Gn.deletevent(failevent)
    Gl=diagnoser(G,failevent,"GL")
    xm=[]
    for state in Gl.X:
        if state[-1]=='Y':
            xm.append(state)
    Gl = Gl.setpar(Xm=xm)
    Gf = simplify(coac(Gl))

    Table,Map = Ri(Gn.Sigma,sigmas[0],1)
    Gni = Gn.renamevents(Map)
    Gni = simplify(Gni.setpar(table=Table))

    #Gni = Gn1//Gn2//Gn3....
    for i in range(1,len(sigmas)):
        Table,Map = Ri(Gn.Sigma,sigmas[i],i+1)
        gni = Gn.renamevents(Map)
        gni = simplify(gni.setpar(table=Table))
        Gni = simplify(Gni//gni)
    G_v = Gni//Gf
    G_v.setgraphic('observer')#Format drawing appearence
    return G_v
