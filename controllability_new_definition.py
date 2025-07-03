from ti_deslab import *

#verificador de controlabilidade

#automato G
syms('a1 a2 b1 b2 1 2 3 4 5 6 7 8')
X = [0, 1, 2, 3, 4, 5, 6, 7, 8]
E =[a1,a2,b1,b2]
sigcontr = [a1,a2,b2]

#tempos em que os eventos sao nao controlaveis
uncontr_dict = { b1 : P.closed(2,5)}

X0 = [0]
Xm = X
T = [(0,a1,1),(0,a2,3),
     (1,b1,2),(1,a2,4),
     (2,a2,5),
     (3,a1,4),(3,b2,6),
     (4,b1,5),(4,b2,7),
     (5,b2,8),
     (6,a1,7),
     (7,b1,8)]

mu = {(0,a1,1) : P.closed(0,2),
      (0,a2,3) : P.closed(2,5),
      (1,b1,2) : P.closed(1,4),
      (1,a2,4) : P.closed(2,4),
      (2,a2,5) : P.closed(1,9),
      (3,a1,4) : P.closed(0,3),
      (3,b2,6) : P.closed(1,7),
      (4,b1,5) : P.closed(0,5),
      (4,b2,7) : P.closed(1,5),
      (5,b2,8) : P.closed(0,5),
      (6,a1,7) : P.closed(1,3),
      (7,b1,8) : P.closed(3,7)
     }

Gtest = fsa(X,E,T,X0,Xm,Sigcon=sigcontr,name="$G$")
G = tia(Gtest,mu)
ti_draw(G)




#automato H
Xh = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
#E =[a1,a2,b1,b2]
#sigcontr = [a1,a2]
X0h = [0]
Xmh = Xh
Th = [(0,a1,1),(0,a2,3),
     (1,b1,2),(1,a2,9),
     (2,a2,5),
     (3,a1,4),(3,b2,6),
     (9,b1,5),(4,b2,7),
     (5,b2,8),
     (6,a1,7),
     (7,b1,8)]

muh = {(0,a1,1) : P.closed(1,2),
       (0,a2,3) : P.closed(2,4),
       (1,b1,2) : P.closed(2,3),
       (1,a2,9) : P.closed(3,4),
       (2,a2,5) : P.closed(1,2),
       (3,a1,4) : P.closed(1,2),
       (3,b2,6) : P.closed(2,3),
       (9,b1,5) : P.closed(0,2),
       (4,b2,7) : P.closed(1,4),
       (5,b2,8) : P.closed(1,4),
       (6,a1,7) : P.closed(2,3),
       (7,b1,8) : P.closed(3,7)
      }

Gtesth = fsa(Xh,E,Th,X0h,Xmh,Sigcon=sigcontr,name="$H$")
H = tia(Gtesth,muh)
ti_draw(H)



#automato Ha = HxG
Hat = ti_product(H,G)
Hat = tia(Hat[0].setpar(name = '$H_a = Prod(H,G)$'),Hat[1])

ti_draw(Hat)

Ha = Hat[0]
Ha_mu = Hat[1]
Ha_mu_new = Hat[1].copy()
NewMarked = list(Hat[0].X) + ['NEW']
Ha = Ha.addstate('NEW')
Ha = Ha.setpar(Xm=NewMarked)



#completar transicoes de Ha
#eestados que chegam eventos uc marcados
contr = G[0].Sigcon
count = 0

for x in Ha.X:
    trans = []
    active_uc = []
    for t in Ha.transitions():
        if (t[0] == x) and (t[1] not in contr):
            trans+=[t]
            active_uc +=[t[1]]

    for ev in active_uc:
            interval = uncontr_dict[ev]
            union = P.empty()
            for t in trans:
                if t[1] == ev:
                    #atualizar intervalo da transicao novo = antigo - intervalo nao contr
                    Ha_mu_new[t] = Ha_mu_new[t] - interval
                    #E SE ESSE INTERVALO DER VAZIO?
                    if len(Ha_mu_new[t]) == 0:
                        Ha = Ha.deletetransition(t)

                    
                    #union = union | Ha_mu_new[t]
            #comp = interval - union
            # no portion, quando subtrai intervalos e o resultado é vazio, nao resulta no P.empty
            # por isso eu checo o comprimento do intervalo
            # o length de qualquer intervalo nao vazio é igual a 1
            # o length de P.empty nao existe
            #if len(comp) != 0 :

            #DÚVIDA: o que fazer com a transicao que acaba ficando com intervalo vazio??
            newt = (x,ev,'NEW')
            count += 1
            Ha_mu_new[newt] = interval
            Ha = Ha.addtransition(newt)
            Ha_new = tia(Ha,Ha_mu_new)
    
    for ev in list(set(G[0].Sigma) - set(active_uc) - set(G[0].Sigcon)):
        if x != 'NEW':
            newt = (x,ev,'NEW')
            Ha_mu_new[newt] = uncontr_dict[ev]
            Ha = Ha.addtransition(newt)
            Ha_new = tia(Ha,Ha_mu_new)
            
if count == 0:
    
    Ha_new = Hat
    Ha_new = tia(Ha_new[0].setpar(name = '$H_{a_{new}}$'),Ha_new[1])
    #Ha_new = tia(Ha_new[0].setpar(Xm=[]),Ha_new[1])

else:  
    Ha_new = tia(Ha_new[0].setpar(name = '$H_{a_{new}}$'),Ha_new[1])
    New_marked = ['NEW']

    #for t in Ha_new[0].transitions():
        #if t[1] not in sigcontr:
            #New_marked += [t[2]]

    #tirar a parte acessivel por causa das transicoes deletadas
    Ha_new = tia(ac(Ha_new[0].setpar(Xm = New_marked)),Ha_new[1])

ti_draw(Ha_new)



Hi = ti_product(Ha_new,G)
Hi = tia(Hi[0].setpar(name = '$H_i = Prod(H_{a_{new}},G)$'),Hi[1])
ti_draw(Hi)


if isitempty(coac(Hi[0])) == True:
    print('É controlável')
else:
    print('Não é controlável')











            
        
    
