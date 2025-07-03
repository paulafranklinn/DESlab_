#    Copyright (C) 2011-2012 by
#    Leonardo Bermeo <lbermeoc@unal.edu.co>
#    Joao Carlos Basilio <basilio@poli.ufrj.br>
#    BSD license.
from deslab import *

def supCont(H,G):
    """
    - This function computes automaton Hi such that Lm(Hi) is the supremal controllable sublanguage
    of Lm(H) with respect to L(G) and Euc = (G.Sigma - G.Sigcon)
    - WARNING: Automaton H must be a nonblocking automaton
    
    """
    Gm = G.setpar(Xm = G.X)
    Euc = G.Sigma - G.Sigcon
    
    Hi = product(H,Gm,simplify=False)
    Hi = Hi.setpar(Sigcon = G.Sigcon, Sigobs = G.Sigobs, name='supC(L(%s))'%(H.name))
    aux = 1
    while aux:
        aux = 0
        for (x,xg) in Hi.X:
            if not(Gm.Gamma(xg) & Euc <= Hi.Gamma((x,xg))):
                Hi = Hi.deletestate((x,xg))
                aux = 1
        Hi = trim(Hi)
        if Hi== fsa():
            return Hi
    Hi = Hi.renamestates('number')
    return Hi

def is_cont(H,G):
    """
    - This function verifies the controllability
    of Lm(H) with respect to L(G) and Euc = (G.Sigma - G.Sigcon)
    - WARNING: Automaton H must be a nonblocking automaton
    
    """
    Gm = G.setpar(Xm = G.X)
    Euc = G.Sigma - G.Sigcon
    
    Hi = product(H,Gm,simplify=False)
    Hi = Hi.setpar(Sigcon = G.Sigcon, Sigobs = G.Sigobs, name='supC(L(%s))'%(H.name))
    for (x,xg) in Hi.X:
        if not(Gm.Gamma(xg) & Euc <= Hi.Gamma((x,xg))):
            return False
    return True
