using JLD2
using FastGaussQuadrature
using HCubature
##

# asymmetric cubature for a complete polynomial system of degree n.
#n = 3
#order = 2*n-1 
#xa, wa = gausslegendre(n) 
#nodes, weights = tensorrule(xa, wa, xa, wa, 2)
nodes, weights = asymmetriccubature(12)
#nodes = [0.9794419976218225 0.9633990795338515; -0.9703748104412847 -0.7627189848126608; 0.6760616605896528 0.01526402220843831; -0.9519046194815316 -0.37216214990235796; -0.9854768681112728 -0.1532257150656241; -0.9566109200707369 0.3009095375359215; 0.9614753684570252 -0.5386257999255943; -0.9712256187052911 0.7267075121245229; -0.9854848813861774 -0.982762018373153; -0.8621183964001193 -0.9288657828332951; -0.9781904628439096 0.9755064530655528; -0.834381386031461 -0.5710992774243813; 0.7251946083497552 0.6249772322191717; 0.9041901730323861 0.8273818040703205; 0.07805423623882504 0.5569334243909445; -0.8115005157616206 0.5660170143572506; 0.7603658826104583 0.9722090864319872; -0.8618359406701791 0.902743794600344; -0.6071932445536927 -0.9929303934469247; -0.6482589438664801 -0.8016144041691429; 0.9111244284300749 -0.9787859099773635; -0.6127234525577426 -0.316686039306279; -0.8261588846535368 -0.011516957931922575; -0.5990375095820828 0.29777391513929374; 0.881944594209698 0.30455423168401186; -0.5943532288812047 0.7993024666081291; -0.6639191806845931 0.990479594832103; -0.33595896050499147 -0.9287099684093596; 0.9860496959216822 -0.8741402356092595; -0.36424565934170633 -0.6100804864228875; 0.9843451674136828 0.5758720461627919; -0.34624194914086687 -0.001282056088579724; 0.8351339001631602 -0.2969215295373109; -0.31970314839220293 0.5890259559968429; 0.5295333676678639 0.8677332339767166; -0.3064575994318151 0.940162544352925; 0.05395502067243376 -0.9706102782272886; -0.06317510961808614 -0.8149202400312969; 0.16083418144353412 -0.5944049250623382; -0.08974561639578102 -0.3118848861738113; 0.9740033846023897 -0.032711819158454024; -0.04465606269381008 0.30654445078980525; 0.6391889604863805 -0.6069346022643242; -0.0012656898874159911 0.8079546049385777; 0.06018192129517237 0.9724620889406334; 0.3435558437204056 -0.9864964778582583; 0.38467295207720137 -0.8334727957413134; 0.8477737460844826 -0.7841024568576497; 0.4376151437462547 -0.3227643212514152; 0.21843220360575388 0.0011012005528191577; 0.49633768035298476 0.333018804515583; 0.2653960225022006 0.5690974851272317; 0.328311898479548 0.7335391252718; 0.29598107126503864 0.9739112095599987; 0.6666628267338494 -0.9427041038894781]
#weights = [0.007678685646670747, 0.026726107101009533, 0.1289853634265132, 0.035919237023755796, 0.018325457153280874, 0.05146780339901072, 0.04039254919407308, 0.027725341927841596, 0.004002012005499107, 0.03470335651613574, 0.006071097399634351, 0.0799608278698593, 0.10679099594239734, 0.04642708616049801, 0.013404666351836786, 0.09257590261672619, 0.026343369758911564, 0.04023637598057319, 0.016623231351755988, 0.08521617508207148, 0.012696232280240539, 0.1455895074546349, 0.11251818503587646, 0.14795895969882147, 0.08946461122197738, 0.09323893188308882, 0.018037880951428373, 0.05665292731460488, 0.012101889252377136, 0.1390739932216903, 0.025400970898628392, 0.1718790344445502, 0.09613149401342268, 0.14818199169586568, 0.07647537308195793, 0.05887632817794223, 0.03166537228575987, 0.08428044175915404, 0.13026311938755222, 0.1693116392364613, 0.039867124971170925, 0.17085089720018198, 0.11104173085396314, 0.10357715416572365, 0.018879239606843446, 0.015498677677690547, 0.09120619103028234, 0.058772143900917596, 0.1465184274100114, 0.17420648907434932, 0.13932040092310008, 0.105748292407316, 0.043344980243863816, 0.028241215740579265, 0.043552509589915654]
order = 12
nodes2, weights2 = contnonsymmetricquad(nodes, weights, order)
##
f(x, y) = x^11*y^4
println(length(weights2))
println(sum([f(nodes2[i, 1], nodes2[i,2])* weights2[i] for i = 1:length(weights2)]))
hcubature(x->f(x[1], x[2]), [-1,-1], [1, 1])

using Plots
scatter(nodes[:, 1], nodes[:,2])

sum(nodes .* weights)
sum(nodes)
