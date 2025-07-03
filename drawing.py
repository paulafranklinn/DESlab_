"""drawing methods"""
#    Copyright (C) 2011-2012 by
#    Leonardo Bermeo <lbermeoc@unal.edu.co>
#    Joao Carlos Basilio <basilio@poli.ufrj.br>
#    GNU license.

#import  deslab
import subprocess
from deslab.src.exceptions import deslabError
from deslab.src.def_const import *
import networkx as nx
import re
from subprocess import call, check_call, check_output, CalledProcessError
import os, sys
import inspect
from deslab.graphics.working.dot2tex_deslab import convert_graph as latex_code
import warnings

VIEWERS             = {'evince':'evince', 'acrobat reader':'acroread'}
VIEWER              = VIEWERS['evince']
DOTINTERFACE        =  'DotInterfaceFile.dot' 
BEAMERPREAMBLE      =  'BeamerPreamble.tex'
PAGEPREAMBLE        =  'PagePreamble.tex'
FIGUREPREAMBLE      =  'FigurePreamble.tex'
TEXPAGEOUT          =  'TexOutput.tex'
WORKING             =  'working'
OUTPUT              =  'output'
TEXFILES            =  'texfiles'
patternDim =re.compile(r'\\node \(\w\d+\) at \((?P<coordX>\d+\.?\d*)pt,(?P<coordY>\d+\.?\d*)pt\)')
dir_path = {WORKING:'', OUTPUT:'',TEXFILES:''}


# CLASS for graphic options

STATE_LAYOUT = {'normal': {'state':'inner sep= 0.25pt, minimum size=0pt, circle,','initpos':''},
    'rectangle': {'state': 'minimum height=0pt, inner sep=0.3pt, inner xsep=0.1pt, rectangle', 'initpos': ''}, 
    'crectangle': {'state': 'minimum height=0mm, inner sep=2mm, chamfered rectangle', 'initpos': ''}, 
    'verifier': {'state': 'minimum height=0pt, inner sep=0.3pt, inner xsep=0.1pt, rectangle', 'initpos': 'above'},  
    'diagnoser': {'state': 'minimum height=0pt, inner xsep=0.1pt, inner ysep=0.3pt, rectangle', 'initpos': 'above'}, 
    'observer': {'state': 'minimum height=0pt, inner sep=0.3pt, inner xsep=0.1pt, rectangle', 'initpos': 'above'}, 
    'vertical': {'state': 'inner sep=0.2pt, minimum size=0pt, circle', 'initpos': 'above'}}  

class graphic:
    
    def __init__(self, style = 'normal', program = 'dot', ranksep = 0.25, nodesep= 0.25, direction = 'LR',
                 FillColor=('plantfill',76), LineColor= ('plantline',85)):
        
        
        self.style = style
        self.program    = program
        self.ranksep    = ranksep
        self.nodesep    = nodesep
        
        try:
            self.state = STATE_LAYOUT[style]['state']
            self.initpos    = STATE_LAYOUT[style]['initpos'] 
        except:
            pass
        
        if style == 'rectangle':
            self.direction  = 'LR'
            self.FillColor  = ('plantfill',76)
            self.LineColor  = ('plantline',85)               
            
        elif style == 'verifier':
            self.direction  = 'LR'
            self.FillColor  = ('yellowfill',76)
            self.LineColor  = ('yellowline',85)      
        
        elif style == 'diagnoser':
            self.direction  = 'UD'
            self.FillColor  = ('skyfill',76)
            self.LineColor  = ('skyline',85)    
        
        elif style == 'observer':
            self.direction  = 'UD'
            self.FillColor  = ('skyfill',76)
            self.LineColor  = ('skyline',85)      
        
        elif style == 'crectangle':
            self.direction  = 'LR'
            self.FillColor  = ('superfill',76)
            self.LineColor  = ('superline',85)    
            
        elif style == 'vertical':
            self.direction  = 'UD'
            self.FillColor  = ('superfill',76)
            self.LineColor  = ('superline',85)       
        
        else:
            self.direction = direction
            self.FillColor = FillColor
            self.LineColor  = LineColor
            self.state = STATE_LAYOUT['normal']['state']
            self.initpos   = STATE_LAYOUT['normal']['initpos']
            
        return
        
 

# TEMPLATES FOR FIGURES

FIGURE_TEMPLATE=r"""\documentclass{beamer}
% basic packages
\usepackage{xcolor}
\usepackage[utf8]{inputenc}

\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{lmodern}
\usepackage[T1]{fontenc}
% tikz libraries
\usepackage{tikz}
\usetikzlibrary{decorations,arrows,shapes,automata,shadows}
\usetikzlibrary{decorations.markings}
% crop preview environment
%\usepackage[active,tightpage]{preview}
%\PreviewEnvironment{tikzpicture}
%\setlength\PreviewBorder{0pt}%
%colors
\definecolor{plantfill}{rgb}{0.960784,   0.850980,   0.039216}
\definecolor{plantline}{rgb}{0.77647,   0.53725,   0.00000}
\definecolor{superfill}{rgb}{0.54510,   0.88235,   0.15686}
\definecolor{superline}{rgb}{0.227451,   0.486275 ,  0.054902}
\definecolor{specfill}{rgb}{0.42353,   0.43922,   0.72157}
\definecolor{specline}{rgb}{0.00000,   0.05098,   0.36471}
\definecolor{uobsb}{rgb}{0.10980,   0.30196,   0.94510 } 
\definecolor{uobsg}{rgb}{0.035294,   0.533333,   0.458824 } 
\definecolor{auto_color}{rgb}{0,0,0}
\definecolor{bg_color}{rgb}{ 0.98823529,  0.98823529,  0.58039216}
\definecolor{skyfill}{rgb}{  0.50196,   0.70196,   1.00000}
\definecolor{skyline}{rgb}{0.00000,   0.23922,   0.75294 }
\definecolor{yellowfill}{rgb}{1.00000,   0.90196,   0.50196}
\definecolor{yellowline}{rgb}{1.00000,   0.4,   0}

% special arrow definition
\newdimen\dima
\newdimen\dimb

\pgfarrowsdeclare{deslab}{deslab}
{
  \dima=0.05pt%
  \advance\dima by.3\pgflinewidth%
  \dimb=6\dima\advance\dimb by.5\pgflinewidth%
  \pgfarrowsleftextend{+-\dimb}
  \dimb=2\dima\advance\dimb by0.5\pgflinewidth%
  \pgfarrowsrightextend{+\dimb}
}
{
  \dima=0.05pt%
  \advance\dima by.3\pgflinewidth%
  \pgfsetdash{}{-0pt}
  \pgfsetroundjoin
  \pgfpathmoveto{\pgfqpoint{2\dima}{0\dima}}
  \pgfpathcurveto
  {\pgfqpoint{-.5\dima}{.5\dima}}
  {\pgfqpoint{-3\dima}{1.5\dima}}
  {\pgfqpoint{-6\dima}{3.25\dima}}
  \pgfpathcurveto
  {\pgfqpoint{-3\dima}{1\dima}}
  {\pgfqpoint{-3\dima}{-1\dima}}
  {\pgfqpoint{-6\dima}{-3.25\dima}}
  \pgfpathcurveto
  {\pgfqpoint{-3\dima}{-1.5\dima}}
  {\pgfqpoint{-.5\dima}{-.5\dima}}
  {\pgfqpoint{2\dima}{0\dima}}
  \pgfpathclose
  \pgfusepathqfillstroke
}

\pgfarrowsdeclare{deslab'}{deslab'}
{
  \dima=0.05pt%
  \advance\dima by.3\pgflinewidth%
  \pgfarrowsleftextend{+-4\dima}
  \pgfarrowsrightextend{+6\dima}
}
{
  \dima=0.05pt%
  \advance\dima by.3\pgflinewidth%
  \pgfpathmoveto{\pgfqpoint{6\dima}{0\dima}}
  \pgfpathcurveto
  {\pgfqpoint{3.5\dima}{.5\dima}}
  {\pgfqpoint{-1\dima}{1.5\dima}}
  {\pgfqpoint{-4\dima}{3.75\dima}}
  \pgfpathcurveto
  {\pgfqpoint{-1.5\dima}{1\dima}}
  {\pgfqpoint{-1.5\dima}{-1\dima}}
  {\pgfqpoint{-4\dima}{-3.75\dima}}
  \pgfpathcurveto
  {\pgfqpoint{-1\dima}{-1.5\dima}}
  {\pgfqpoint{3.5\dima}{-.5\dima}}
  {\pgfqpoint{6\dima}{0\dima}}
  \pgfusepathqfill
}

\pgfarrowsdeclare{init}{init}
{
  \dima=0.05pt%
  \advance\dima by.25\pgflinewidth%
  \dimb=5.5\dima\advance\dimb by.5\pgflinewidth%
  \pgfarrowsleftextend{+-\dimb}
  \dimb=.5\dima\advance\dimb by0.707\pgflinewidth%
  \pgfarrowsrightextend{+\dimb}
}
{
  \dima=0.05pt%
  \advance\dima by.25\pgflinewidth%
  \pgfsetdash{}{+0pt}
  \pgfsetmiterjoin
  \pgfpathmoveto{\pgfqpoint{-5.5\dima}{-6\dima}}
  \pgfpathlineto{\pgfqpoint{0.5\dima}{0\dima}}
  \pgfpathlineto{\pgfqpoint{-5.5\dima}{6\dima}}
  \pgfpathclose
  \pgfusepathqfillstroke
}

\pgfarrowsdeclare{init'}{init'}
{
  \dima=0.05pt%
  \advance\dima by.25\pgflinewidth%
  \pgfarrowsleftextend{+-.5\pgflinewidth}
  \dimb=6\dima\advance\dimb by0.707\pgflinewidth%
  \pgfarrowsrightextend{+\dimb}
}
{
  \dima=0.05pt%
  \advance\dima by.25\pgflinewidth%
  \pgfsetdash{}{+0pt}
  \pgfsetmiterjoin
  \pgfpathmoveto{\pgfqpoint{0\dima}{-6\dima}}
  \pgfpathlineto{\pgfqpoint{6\dima}{0\dima}}
  \pgfpathlineto{\pgfqpoint{0\dima}{6\dima}}
  \pgfpathclose
  \pgfusepathqstroke
}
% start here
\begin{document}
% setting preview
\pagestyle{empty}
\enlargethispage{100cm}
% font selection
\fontsize{1pt}{1pt}
\selectfont
\tikzset{obs_edge arrow/.style={->, >=deslab, line width =0.075pt}}
\tikzset{accepting/.style={double distance=0.1pt}}

"""

BEAMER_TEMPLATE = r"""\documentclass{beamer}
% basic packages
\usepackage{xcolor}
\usepackage[utf8]{inputenc}

\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{lmodern}
\usepackage[T1]{fontenc}
% tikz libraries
\usepackage{tikz}
\usetikzlibrary{decorations,arrows,shapes,automata,shadows}
\usetikzlibrary{decorations.markings}
%colors
\definecolor{plantfill}{rgb}{0.960784,   0.850980,   0.039216}
\definecolor{plantline}{rgb}{0.77647,   0.53725,   0.00000}
\definecolor{superfill}{rgb}{0.54510,   0.88235,   0.15686}
\definecolor{superline}{rgb}{0.227451,   0.486275 ,  0.054902}
\definecolor{specfill}{rgb}{0.42353,   0.43922,   0.72157}
\definecolor{specline}{rgb}{0.00000,   0.05098,   0.36471}
\definecolor{uobsb}{rgb}{0.10980,   0.30196,   0.94510 } 
\definecolor{uobsg}{rgb}{0.035294,   0.533333,   0.458824 } 
\definecolor{auto_color}{rgb}{0,0,0}
\definecolor{bg_color}{rgb}{ 0.98823529,  0.98823529,  0.58039216}
\definecolor{skyfill}{rgb}{  0.50196,   0.70196,   1.00000}
\definecolor{skyline}{rgb}{0.00000,   0.23922,   0.75294 }
\definecolor{yellowfill}{rgb}{1.00000,   0.90196,   0.50196}
\definecolor{yellowline}{rgb}{1.00000,   0.4,   0}

% beamer data
\author{LCA--Lab. de Controle e Automac\~ao}
\title{DESLab Software Package}
\usetheme{Madrid}
\usecolortheme{seagull}
\setbeamercolor{normal text}{bg=bg_color!30}
\usefonttheme{professionalfonts}

% special arrow definition
\newdimen\dima
\newdimen\dimb

\pgfarrowsdeclare{deslab}{deslab}
{
  \dima=0.05pt%
  \advance\dima by.3\pgflinewidth%
  \dimb=6\dima\advance\dimb by.5\pgflinewidth%
  \pgfarrowsleftextend{+-\dimb}
  \dimb=2\dima\advance\dimb by0.5\pgflinewidth%
  \pgfarrowsrightextend{+\dimb}
}
{
  \dima=0.05pt%
  \advance\dima by.3\pgflinewidth%
  \pgfsetdash{}{-0pt}
  \pgfsetroundjoin
  \pgfpathmoveto{\pgfqpoint{2\dima}{0\dima}}
  \pgfpathcurveto
  {\pgfqpoint{-.5\dima}{.5\dima}}
  {\pgfqpoint{-3\dima}{1.5\dima}}
  {\pgfqpoint{-6\dima}{3.25\dima}}
  \pgfpathcurveto
  {\pgfqpoint{-3\dima}{1\dima}}
  {\pgfqpoint{-3\dima}{-1\dima}}
  {\pgfqpoint{-6\dima}{-3.25\dima}}
  \pgfpathcurveto
  {\pgfqpoint{-3\dima}{-1.5\dima}}
  {\pgfqpoint{-.5\dima}{-.5\dima}}
  {\pgfqpoint{2\dima}{0\dima}}
  \pgfpathclose
  \pgfusepathqfillstroke
}

\pgfarrowsdeclare{deslab'}{deslab'}
{
  \dima=0.05pt%
  \advance\dima by.3\pgflinewidth%
  \pgfarrowsleftextend{+-4\dima}
  \pgfarrowsrightextend{+6\dima}
}
{
  \dima=0.05pt%
  \advance\dima by.3\pgflinewidth%
  \pgfpathmoveto{\pgfqpoint{6\dima}{0\dima}}
  \pgfpathcurveto
  {\pgfqpoint{3.5\dima}{.5\dima}}
  {\pgfqpoint{-1\dima}{1.5\dima}}
  {\pgfqpoint{-4\dima}{3.75\dima}}
  \pgfpathcurveto
  {\pgfqpoint{-1.5\dima}{1\dima}}
  {\pgfqpoint{-1.5\dima}{-1\dima}}
  {\pgfqpoint{-4\dima}{-3.75\dima}}
  \pgfpathcurveto
  {\pgfqpoint{-1\dima}{-1.5\dima}}
  {\pgfqpoint{3.5\dima}{-.5\dima}}
  {\pgfqpoint{6\dima}{0\dima}}
  \pgfusepathqfill
}

\pgfarrowsdeclare{init}{init}
{
  \dima=0.05pt%
  \advance\dima by.25\pgflinewidth%
  \dimb=5.5\dima\advance\dimb by.5\pgflinewidth%
  \pgfarrowsleftextend{+-\dimb}
  \dimb=.5\dima\advance\dimb by0.707\pgflinewidth%
  \pgfarrowsrightextend{+\dimb}
}
{
  \dima=0.05pt%
  \advance\dima by.25\pgflinewidth%
  \pgfsetdash{}{+0pt}
  \pgfsetmiterjoin
  \pgfpathmoveto{\pgfqpoint{-5.5\dima}{-6\dima}}
  \pgfpathlineto{\pgfqpoint{0.5\dima}{0\dima}}
  \pgfpathlineto{\pgfqpoint{-5.5\dima}{6\dima}}
  \pgfpathclose
  \pgfusepathqfillstroke
}

\pgfarrowsdeclare{init'}{init'}
{
  \dima=0.05pt%
  \advance\dima by.25\pgflinewidth%
  \pgfarrowsleftextend{+-.5\pgflinewidth}
  \dimb=6\dima\advance\dimb by0.707\pgflinewidth%
  \pgfarrowsrightextend{+\dimb}
}
{
  \dima=0.05pt%
  \advance\dima by.25\pgflinewidth%
  \pgfsetdash{}{+0pt}
  \pgfsetmiterjoin
  \pgfpathmoveto{\pgfqpoint{0\dima}{-6\dima}}
  \pgfpathlineto{\pgfqpoint{6\dima}{0\dima}}
  \pgfpathlineto{\pgfqpoint{0\dima}{6\dima}}
  \pgfpathclose
  \pgfusepathqstroke
}



\begin{document}
%font selection
\fontsize{1pt}{1pt}
\selectfont
\tikzset{obs_edge arrow/.style={->, >=deslab, line width =0.075pt}}
\tikzset{accepting/.style={double distance=0.1pt}}

"""

EMPTY_AUTOMATON = r"""\tikzset{every initial by arrow/.style=
{>=initarrow,initial distance=50pt, draw=specline, line width=4pt, initial where= above, initial text={}}}
\fontsize{30pt}{30pt}
\selectfont
\node (s1) at (2,2) [ellipse,line width=1pt, ball color=superline, path fading=north, fill=superfill!70 ] {
\begin{tabular}{c}
\\
{\textsc{Empty }}
{\textsc{Automaton}}\\
\hline
\\
\textcolor{white}{ Automaton without states}\\
\\
\end{tabular}};

"""

# DESIGN OF STATE RENDERING 

#STATE_LAYOUT = {'normal': {'state':'inner sep= 0.25pt, minimum size=0pt, circle,','init_pos':''},
#                'rectangle': {'state': 'minimum height=0pt, inner sep=0.3pt, inner xsep=0.1pt, rectangle', 'init_pos': ''}, 
#                'crectangle': {'state': 'minimum height=0mm, inner sep=2mm, chamfered rectangle', 'init_pos': ''}, 
#                'verifier': {'state': 'minimum height=0pt, inner sep=0.3pt, inner xsep=0.1pt, rectangle', 'init_pos': 'above'},  
#                'diagnoser': {'state': 'minimum height=0pt, inner xsep=0.1pt, inner ysep=0.3pt, rectangle', 'init_pos': 'above'}, 
#                'observer': {'state': 'minimum height=0pt, inner sep=0.3pt, inner xsep=0.1pt, rectangle', 'init_pos': 'above'}, 
#                'vertical': {'state': 'inner sep=0.2pt, minimum size=0pt, circle', 'init_pos': 'above'}}   

COLOR_TEMPLATE = ''
PREAMBLE_DIC = {'beamer': BEAMER_TEMPLATE,
                'figure': FIGURE_TEMPLATE,
                'figurecolor': FIGURE_TEMPLATE}

global fig_counter, window_counter 
fig_counter = 1
window_counter = 0



def create_digraph(self) :
    """ This function creates a digraph containing a map
     between an alias for each state, and event, in order
     to draw a right graph with standarized label size  """
     
    def tex(symbol):
        """ This function converts the symbol used in textual
        interface to the latex label that will be renderized
        """
        if self.symDict == {}:
            converted = symbol        
        elif symbol in self.symDict:
            converted = self.symDict[symbol]
        elif symbol in self.X | self.Sigma:
            converted = symbol 
        else: 
            raise invalidLabel('Symbol %s was not recognized'%(symbol))
        return converted  
           
    Graph       = self.Graph    
    nxDigraph   = nx.MultiDiGraph()
    Siguobs     = self.Sigma-self.Sigobs
    Sigobs      = self.Sigobs
    Sigcon      = self.Sigcon  
    Siguncon    = self.Sigma-self.Sigcon 
    
    if os.name == 'nt':
        styleNodeDic ={'initial': '"state,initial"',
                       'marked': '"state,accepting"',
                       'initialmarked': '"state,initial,accepting"',
                       'normal':'"state"'} 
    else: 
        styleNodeDic ={'initial': "state,initial",
                       'marked': "state,accepting",
                       'initialmarked': "state,initial,accepting",
                       'normal': "state"}       
    for node in Graph.nodes(data=True):
        symbol = node[0]
        if symbol in (self.Xm & self.X0) :
            actualStyle=styleNodeDic['initialmarked']
        elif symbol in self.Xm :
            actualStyle=styleNodeDic['marked']
        elif symbol in self.X0: 
            actualStyle=styleNodeDic['initial']            
        else:
            actualStyle=styleNodeDic['normal']
        nxDigraph.add_node(node[1]['label'], label = tex(symbol), style=actualStyle)
        #nxDigraph.add_node(node[1]['label'], label = tex(symbol), style=actualStyle)
        
    for source in Graph.nodes():               
        for target in Graph.successors(source):
            events       = Graph.get_edge_data(source,target).keys()
            events_obs   = [tex(i)+',' for i in sorted(list(set(events) & Sigobs))]                  
            events_unobs = [tex(i)+',' for i in sorted(list(set(events) & Siguobs))]              
            events_obs   = ''.join(events_obs).rstrip(',')         
            events_unobs = ''.join(events_unobs).rstrip(',')                  
            source_dig   = Graph._node[source]['label']
            target_dig   = Graph._node[target]['label']
            
            if  events_obs != '':                       
                if  os.name == 'nt':   # in the case if windows system               
                    events_obs = r'"'+events_obs+r'"'
                nxDigraph.add_edge(source_dig, target_dig, key=events_obs,
                                   label=events_obs, style='obs_edge arrow') # we have to fix the labels problem
            if  events_unobs != '':               
                if  os.name == 'nt':             
                    events_unobs = r'"'+events_unobs+r'"'  
                nxDigraph.add_edge(source_dig, target_dig, key=events_unobs,
                                   label=events_unobs, style='unobs_edge arrow')
    
   
    return nxDigraph



##def create_digraph_dot(self) : #parece que nao usa
##    """ This function creates a digraph containing a map
##     between an alias for each state, and event, in order
##     to draw a right graph with standarized label size  """
##     
##    def tex(symbol):
##        """ This function converts the symbol used in textual
##        interface to the latex label that will be renderized
##        """
##        if self.symDict == {}:
##            converted = symbol        
##        elif symbol in self.symDict:
##            converted = self.symDict[symbol]
##        elif symbol in self.X | self.Sigma:
##            converted = symbol 
##        else: 
##            raise invalidLabel('Symbol %s was not recognized'%(symbol))
##        return converted  
##           
##    Graph       = self.Graph    
##    nxDigraph   = nx.MultiDiGraph()
##    Siguobs     = self.Sigma-self.Sigobs
##    Sigobs      = self.Sigobs
##    Sigcon      = self.Sigcon  
##    Siguncon    = self.Sigma-self.Sigcon 
##    
##    if os.name == 'nt':
##        styleNodeDic ={'initial': '"state,initial"',
##                       'marked': '"state,accepting"',
##                       'initialmarked': '"state,initial,accepting"',
##                       'normal':'"state"'} 
##    else: 
##        styleNodeDic ={'initial': "state,initial",
##                       'marked': "state,accepting",
##                       'initialmarked': "state,initial,accepting",
##                       'normal': "state"}       
##        
##    for node in Graph.nodes_iter(data=True): 
##        symbol = node[0]       
##        if symbol in (self.Xm & self.X0) :
##            actualStyle=styleNodeDic['initialmarked']
##        elif symbol in self.Xm :
##            actualStyle=styleNodeDic['marked']
##        elif symbol in self.X0: 
##            actualStyle=styleNodeDic['initial']            
##        else:
##            actualStyle=styleNodeDic['normal']                                      
##        nxDigraph.add_node(node[1]['label'], label = tex(symbol), style=actualStyle)
##    
##        
##    for source in Graph.nodes_iter():               
##        for target in Graph.successors_iter(source):
##            events       = Graph.get_edge_data(source,target).keys()
##            events_obs   = [tex(i)+',' for i in sorted(list(set(events) & Sigobs))]                  
##            events_unobs = [tex(i)+',' for i in sorted(list(set(events) & Siguobs))]              
##            events_obs   = ''.join(events_obs).rstrip(',')         
##            events_unobs = ''.join(events_unobs).rstrip(',')                  
##            source_dig   = Graph.node[source]['label']
##            target_dig   = Graph.node[target]['label']
##            
##            if  events_obs != '':                       
##                if  os.name == 'nt':   # in the case if windows system               
##                    events_obs = r'"'+events_obs+r'"'
##                nxDigraph.add_edge(source_dig, target_dig, key=events_obs,
##                                   label=events_obs, style='obs_edge arrow') # we have to fix the labels problem
##            if  events_unobs != '':               
##                if  os.name == 'nt':             
##                    events_unobs = r'"'+events_unobs+r'"'  
##                nxDigraph.add_edge(source_dig, target_dig, key=events_unobs,
##                                   label=events_unobs, style='unobs_edge arrow')
##    
##   
##    return nxDigraph




def setupdir():    
    """This function recover the absolute path of the user and 
    set the basic directories for drawing"""
    path_drawing       = inspect.currentframe().f_code.co_filename
    graphics_path      = os.path.dirname(path_drawing)    
    dir_path[WORKING]  = os.path.join(graphics_path, WORKING)
    dir_path[OUTPUT]   = os.path.join(graphics_path, OUTPUT)
    dir_path[TEXFILES] = os.path.join(graphics_path, TEXFILES)    
    return

def adjust_label_dot(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Corrige os rótulos sem aspas
    with open(filename, 'w') as file:
        for line in lines:
            if "label=" in line and not "label=\"" in line:
                line = line.replace("label=", "label=\"")
                line = line.replace(", style=", "\", style=")
            file.write(line)

def auto2dot(automaton):    
    """ 
    This function creates a dot description
    of the input automaton
    """
    digraph_aut = create_digraph(automaton)
    name_dotfile = os.path.join(dir_path[WORKING], DOTINTERFACE)     
    #nx.write_dot(digraph_aut, name_dotfile)
    nx.drawing.nx_pydot.write_dot(digraph_aut, name_dotfile)
    file_obj  = open(name_dotfile)
    dot_string = file_obj.read()
    dot_string = dot_string.replace('strict','')
    dot_string = dot_string.lstrip()
    file_obj.close()
    return dot_string


def determine_size(texfile, nodes):
    """This function parses the created text file in order 
    to get the ratio for rendering the graphics in Beamer mode"""
    if nodes == 0:
        size='50mm'
        return size
    coordX=[]
    coordY=[]
    
    for match in patternDim.finditer(texfile):
        coordX.append(float(match.group('coordX')))
        coordY.append(float(match.group('coordY')))
    FigureWidth = max(coordX) - min(coordX)
    FigureHeight = max(coordY) - min(coordY)
    
    miliWidth, miliHeight =  FigureWidth*25.4/72.0, FigureHeight*25.4/72.0 
    
    if (miliWidth == 0) & (miliHeight != 0):
        factorSize = 20/miliHeight
        newSize = miliHeight*factorSize
        size = str(newSize)+'mm'    
  
    elif (miliWidth != 0) & (miliHeight == 0):
        size = '100mm' 
    elif (miliWidth == 0) & (miliHeight == 0):
        size = '20mm'
    else:        
        aspectRelation = miliWidth/miliHeight
        if aspectRelation >= 1.57:            
            size = '100mm'
        else:           
           factorSize = 60/miliHeight
           newSize = miliWidth*factorSize
           size = str(newSize)+'mm'    
    return size




def automaton2tikfig(automaton):
        
               
    #environment variables
    dot_init = auto2dot(automaton)  
    graph = automaton.Graph
    graphic = automaton.graphic
    nodesep = graphic.nodesep
    ranksep = graphic.ranksep 
    program = graphic.program
    direction = graphic.direction
    
    preamble = 'rankdir=%s;\nnodesep=%s;\nranksep=%s;\n'%(direction, str(nodesep), str(ranksep))
    for x_i in automaton.X0 :
        preamble += '\t'+graph._node[x_i]['label']+' '+'[style="state"];\n'
    for x_m in automaton.Xm :
        preamble += '\t'+graph._node[x_m]['label']+' '+'[style="state,accepting"];\n'     
    auto_dotfile = dot_init[0:11] + preamble + dot_init[11:]# texcommands
    file = os.path.join(dir_path[WORKING], DOTINTERFACE)      
    fileObj = open(file, 'w')
    fileObj.write(auto_dotfile)
    fileObj.close()
    adjust_label_dot(file)
    command = '%s -Txdot %s | '%(program, DOTINTERFACE) + 'python dot2tex_deslab.py -ftikz --codeonly --texmode math'
    try:
        fig_texcode = check_output(command,shell=True, cwd = dir_path[WORKING], stderr=subprocess.PIPE)
    except CalledProcessError:
        raise deslabError('I could not create a tex file of the automaton object')
        
    return fig_texcode

    
def tex2pdf(tex_filename) :    
    global window_counter
    window_counter += 1    
    pdf_outputname = 'Figure-'+str(window_counter)+str(fig_counter)
    command = '/Library/TeX/texbin/pdflatex' + ' -interaction=batchmode -no-shell-escape'
    command += ' -output-directory ' + dir_path[OUTPUT] + ' -jobname ' + pdf_outputname + ' ' + tex_filename

    try:
        retcode = check_call(command,shell=True, cwd = dir_path[WORKING])  
    except CalledProcessError:
        raise deslabError("There is a trouble with the generation of the pdf file. Check pdflatex installation")
    return pdf_outputname+'.pdf'
     
    
   
def automaton2page(automaton, style):     
    from deslab.src.automatadefs import fsa   
   
    graphic = automaton.graphic
        
    strLineColor = str(graphic.LineColor[0]) + '!' + str(graphic.LineColor[1])
    strFillColor = str(graphic.FillColor[0]) + '!' + str(graphic.FillColor[1])
    state_ly = graphic.state 
    initpos = graphic.initpos 
    print("generating latex code of automaton")
    if automaton == fsa():
        tikz_code = EMPTY_AUTOMATON          
    else:
        tikz_code = "".join(map(chr,automaton2tikfig(automaton)))
        tikz_code = tikz_code + '\n'  
    
    nodes=len(automaton.X)       
      
              
    if style == 'beamer':
        size = determine_size(tikz_code, nodes) 
        init_tex = r"""        
\tikzset{unobs_edge arrow/.style={->, >=deslab, dash pattern=on 0.1pt off 0.15pt, draw = %s,
line width =0.075pt}}

\tikzset{every initial by arrow/.style=
{>=init, initial distance=1.7pt,line width=0.125pt, initial where=%s, fill =%s, draw =%s, initial text={}}}

\tikzset{every state/.style={draw= %s, fill= %s, line width = 0.1, %s}}
        
\begin{frame}       
\frametitle{\textcolor{auto_color}{\textcolor{white}{\framebox{Fig. \; %s.%s }\quad} %s}}
\begin{center}
\resizebox{%s}{!}{
\begin{tikzpicture}
""" % ('uobsg', initpos, 'yellow!60', 'green!30!black',
         strLineColor, strFillColor,   state_ly, str(window_counter),
         str(fig_counter), automaton.name, size)
       
        figure_texcode = init_tex + tikz_code +'\\end{tikzpicture}} \n \\end{center} \n \\end{frame}'
        
        
    elif style == 'figurecolor':
        size = determine_size(tikz_code, nodes)
        init_tex = r"""        
\tikzset{unobs_edge arrow/.style={->,>=deslab, dash pattern=on 0.1pt off 0.15pt, draw = %s,
line width =0.075pt}}

\tikzset{every initial by arrow/.style=
{>=init, initial distance=1.7pt,line width=0.1pt, initial where=%s, fill =%s, draw =%s, initial text={}}}

\tikzset{every state/.style={draw= %s, fill= %s, line width = 0.1, %s}}

\begin{frame}
\begin{center}
\resizebox{%s}{!}{
\begin{tikzpicture}
""" % ('uobsg',  initpos, 'orange!40', 'red!60!black',
         strLineColor, strFillColor,  state_ly, size)
        figure_texcode = init_tex + tikz_code + '\\end{tikzpicture}} \n \\end{center} \n \\end{frame}'
        
        
    elif style == 'figure':
        size = determine_size(tikz_code, nodes)
        init_tex = r"""        
\tikzset{unobs_edge arrow/.style={->, >=deslab, dash pattern=on 0.1pt off 0.15pt, draw = %s,
line width =0.075pt }}

\tikzset{every initial by arrow/.style=
{>=deslab, initial distance=1.7pt,line width=0.075pt, initial where=%s, fill =%s, draw =%s, initial text={}}}

\tikzset{every state/.style={draw= %s, fill= %s, line width = 0.1, %s}}

\begin{frame}
\begin{center}
\resizebox{%s}{!}{
\begin{tikzpicture}
""" % ('black', initpos, 'black', 'black',
         'black', 'white',   state_ly, size)
        figure_texcode = init_tex + tikz_code + '\\end{tikzpicture}} \n \\end{center} \n \\end{frame}'
                
    else :
        raise deslabError('style %s is not defined yet'%style)
    return figure_texcode 


def write_texfile(TexString,TexfileOut):
    TexComplete = TexString+r'\end{document}'
    Texpresfile = os.path.join(dir_path[WORKING], TexfileOut)      
    fileObj = open(Texpresfile, 'w')
    fileObj.write(TexComplete)
    fileObj.close()
    return TexfileOut

def write_texfile_here(TexString):
    TexComplete = TexString+r'\end{document}'
    tex_outputname = 'Figure'+'.tex'
    Texpresfile = os.path.join(os.getcwd(),tex_outputname)      
    fileObj = open(Texpresfile, 'w')
    fileObj.write(TexComplete)
    fileObj.close()

def openviewer(filepdf):
    filepath = os.path.join(dir_path[OUTPUT], filepdf)
    if sys.platform.startswith('darwin'):  # macOS
        subprocess.call(['open', filepath])
    elif os.name == 'nt':  # Windows
        os.startfile(filepath)
    else:  # Linux/Unix
        subprocess.call([VIEWER, filepath])
    return  
      
def draw(*automatavars, save_file=False):      
    from deslab.src.automatadefs import fsa
    global fig_counter 
    setupdir() # we setup the directories in the current operative system 
    
    # the last argument is an fsa then we use the default style
    if isinstance(automatavars[-1], fsa):          
        AutomataList=list(automatavars[0:-1])
        AutomataList.append(automatavars[-1])
        style='beamer'
    # the last argument  is a list then we have a vector of automata
    # we must setup the default style
    elif isinstance(automatavars[-1],list):   
        AutomataList=automatavars[-1] 
        style='beamer'
    # the last argument  is a style
    # and we have a list of automata   
    elif isinstance(automatavars[-1],str) & isinstance(automatavars[-2],list) :
        AutomataList = automatavars[-2]
        style = automatavars[-1]  
    # we have a comma separated automata and one defined style
    else :
        style=automatavars[-1] 
        AutomataList=list(automatavars[0:-1])
    # we setup the preamble
    if style in PREAMBLE_DIC:
        preamble_tex = PREAMBLE_DIC[style]
    else:
        preamble_tex = PREAMBLE_DIC['beamer']
        style = 'beamer'
   
    fig_counter = 0 
    # we generate the latex code for every automata passed to draw function 
    for automaton in AutomataList: 
        #in the case of a large automata we avoid rendering
        if len(automaton.X)>100:
            print("Automaton ",automaton.name," : Current version of DESlab supports drawings up to 100 states.\n Try the method tmx of the automaton or the function save to get the transition matrix of the automaton.")
            #continue
        else:
            preamble_tex += automaton2page(automaton, style)
        fig_counter += 1
    write_texfile(preamble_tex,TEXPAGEOUT)
    if save_file:
        write_texfile_here(preamble_tex)
    page_pdf = tex2pdf(TEXPAGEOUT)    
    openviewer(page_pdf)
    return 







