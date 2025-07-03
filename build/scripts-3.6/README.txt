###  English ###

========================================================================================
Before installing DESlab 
========================================================================================

Make sure to verify the following steps:
	
	- Install a LaTeX distribution (TeXlive, MikTeX or another). The TeXlive distribution (included in folder "programas" of DESlab package or available at https://www.tug.org/texlive/acquire.html) has been fully tested for this version of DESlab, but other distributions should also work.
	
	- Add the path of pdflatex.exe of TexLive (usually, C:\...\Texlive\2017\bin\win32) to the Windows Environment Variables. (Tutorial below)

	- Install o python 3.6 (included in folder "programas" of DESlab package or available at https://www.python.org/ftp/python/3.6.4/python-3.6.4-amd64.exe).
	  Remark: During the installation of python 3.6, there should be a checkbox to include its path in the Windows Environment Variables.
	  
	- If you have an old version of python installed, remove its path from the Windows Environment Variables, and add the path of python 3.6. (Tutorial below)


========================================================================================
Setting the Windows Environment Variables
========================================================================================

	- Select System from the Control Panel, selecting Advanced system settings, and clicking Environment Variables.
	
	- In the Environment Variable window, edit variable PATH in the System Variables. In  the variable value, the several paths are joined by using operator ';' (semicolon). Delete the path of the older version of python and insert the path of python 3.6.
	Remark: the path should be path to the folder, not including the file, for instance, 
		- WRONG -> C:/...../Python/Pytohn36-32/python.exe;
		- RIGHT -> C:/...../Python/Pytohn36-32;


========================================================================================
DESlab Installation
========================================================================================

1) Execute script "Install.bat" at folder "DESlab". The following modules and software are installed:
	- Graphviz 2.28
	- FaDo 1.3.5.1
	- Future 0.16.0
	- Networkx 2.1
	- Pyparsing 2.2.0
	- Pydot 1.2.4	
	- DESlab
	
2) Verify if the path of Graphviz (usually, ...Graphviz\bin) has been added to the Windows Environment Variables.

3) Run the following commands in python IDLE to test the installation:

from  deslab import *
syms('q1 q2 q3 a1 b1 e f')
table = [(a1,'a_1'),(b1,'b_1'),(q1,'q_1'),(q2,'q_2'),(q3,'q_3')]
X = [q1,q2,q3]
Sigma = [a1,b1,e]
X0 = [q1]
Xm = [q3]
T =[(q1,b1,q2),(q2,b1,q3),(q3,e,q3)]
G1 = fsa(X,Sigma,T,X0,Xm,table,name='$G_1$')
draw(G1)
draw(G1,'figurecolor')
