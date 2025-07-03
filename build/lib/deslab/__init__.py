"""
DESlab
========

    DESlab a scientific computing program written in python
    for the development of algorithms for analysis and synthesis
    of discrete event systems (DES) modeled as automata.
    The main objective of DESlab is to provide a unified tool 
    that integrates automata, graph algorithms, and numerical
    calculations. DESlab also allows the definition of symbolic
    variables of type automaton and incorporates concise
    instructions to manipulate, operate, analyze and visualize
    these variables, with a syntax and an abstraction level close
    to the notation used in DES theory

    http://www.dee.ufrj.br/controle_automatico/

Using
-----

    Just write in Python

    >>> from deslab import *
    >>> syms('q1 q2 q3 a1 b1 e')
    >>> table = [(a1,'\alpha_1'),(b1,'\beta_1'),(q1,'q_1'),(q2,'q_2'),(q3,'q_3')]
    >>> X = [q1,q2,q3]
    >>> Sigma = [a1,b1,e]
    >>> X0 = [q1]
    >>> Xm = [q1,q3]
    >>> T =[(q1,b1,q2),(q2,b1,q3),(q2,e,q3),(q3,a1,q1),(q2,a1,q2)]
    >>> G1 = fsa(X,Sigma,T,X0,Xm,table,name='$G_1$ -- example 1')
    >>> draw(G1)
"""
#    Copyright (C) 2011-2012 by
#    Leonardo Bermeo <lbermeoc@unal.edu.co>
#    Joao Carlos Basilio <basilio@poli.ufrj.br>
#    BSD license.


# adding source code
import deslab.src
from deslab.src.automatadefs import *
from deslab.src.utilities import *
from deslab.src.exceptions import *
from deslab.src.algorithms import *
from deslab.src.structure import *
from deslab.src.comparison import *
from deslab.src.graphs import *

# adding graphics
import deslab.graphics
from deslab.graphics.drawing import *

import deslab.graphics.working 
from deslab.graphics.working import *

# adding toolboxes
import deslab.toolboxes
from deslab.toolboxes.diagnosis import *
from deslab.toolboxes.supervisory import *
#from deslab.toolboxes.supervisory import *
from deslab.toolboxes.opacity_verifier import *
from deslab.toolboxes.opacity_enforcement import *
from deslab.toolboxes.ti_functions import *
from deslab.toolboxes.ti_diagnosis import *
from deslab.toolboxes.ti_opacity_verifier import *

# adding readwrite
import deslab.readwrite
from deslab.readwrite.inputoutput import *
