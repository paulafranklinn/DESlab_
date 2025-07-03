#    Copyright (C) 2011-2012 by
#    Leonardo Bermeo <lbermeoc@unal.edu.co>
#    Joao Carlos Basilio <basilio@poli.ufrj.br>
#    GNU license.
import deslab
from deslab.src.def_const import *

          
 


##### These are the classes defining the software and user common errors


class desError(Exception):
    #Base class for exceptions in this module.
    pass

class inputError(Exception):
    #Base class for exceptions in this module.
    pass



class DFAerror(desError):
    pass


class notDFAError(DFAerror):  
        
    def __init__(self,msg):
        self.message =  msg
    def __str__(self):
        return "It is not a DFA. %s"%(self.message)

class stateError(desError):   
          
    def __init__(self,msg):
        self.message =  msg
        
    def __str__(self):
        return "Invalid State. %s"%(self.message)

class deslabError(desError):   
          
    def __init__(self,msg):
        self.message =  msg
        
    def __str__(self):
        return self.message
    
class stateMembershipError(desError):   
          
    def __init__(self,msg):
        self.message =  msg
        
    def __str__(self):
        return "Undefined state. %s"%(self.message)

class eventMembershipError(desError):   
          
    def __init__(self,msg):
        self.message =  msg
        
    def __str__(self):
        return "Undefined event. %s"%(self.message)
    
class markedSetError(desError):   
          
    def __init__(self,msg):
        self.message =  msg
        
    def __str__(self):
        return "Error in marked set Xm. %s"%(self.message)
    
class initialStateError(desError):   
          
    def __init__(self,msg):
        self.message =  msg        
    def __str__(self):
        return "Bad definition of X0 state. %s"%(self.message)
    
class epsilonDFAError(DFAerror):  
        
    def __init__(self,msg):
        self.message =  msg
    def __str__(self):
        return "It is not a DFA because it has epsilon events. %s"%(self.message)
    
class inputStringError(desError):  
        
    def __init__(self,msg):
        self.message =  msg
    def __str__(self):
        return "It is a valid input string. %s"%(self.message)
    
class invalidAutomaton(desError):  
        
    def __init__(self,msg):
        self.message =  msg
    def __str__(self):
        return "It is a nont valid automaton. %s"%(self.message)
    
    
class invalidArgument(inputError):
    def __init__(self,msg):
        self.message =  msg
    def __str__(self):
        return "Invalid input argument. %s"%(self.message)
    
class invalidLabel(inputError):
    def __init__(self,msg):
        self.message =  msg
    def __str__(self):
        return "Invalid input argument. %s"%(self.message)
    
class invalidTransition(inputError):
    def __init__(self,msg):
        self.message =  msg
    def __str__(self):
        return "Invalid input argument. %s"%(self.message)
