# Modeling the spread of Xylella Fastidiosa in Apulia, Italy

In this project we adopt the kernel spread model proposed by White et al. (2017, [[1]](#1)) in their study "Modelling the spread and control of Xylella fastidiosa in the early stages of invasion in Apulia, Italy" to estimate the speed of disease spread and to analyze the effect of the implemented emergency measures. 
Through this analysis, we aim to determine whether more effective strategies could have been employed in managing the outbreak.



## What is Xylella?

**_Xylella fastidiosa_** is a bacterial pathogen that is transmitted by insect vectors feeding on xylem sap from trees. 
Xylem cells act as tubes that transport water and nutrients (sap) from the soil up to the branches and leaves. 
The bacterium acts by obstructing the flow of nutrients through the xylem network, causing the tree to slowly die out, starting from the leaves and proceeding up to the main trunk.
Common symptoms are leaf scorching, stunting and wilting, and a reduction in the size and quality of fruit production.
The figures below (taken from [[2]](#2)) capture the symptoms that detail the progress of the disease in the leaves (A), branches (B) and tree canopy (C) in the olive tree and oleander, respectively.  

<div align="center">
    <img src=Images/Symptoms_on_olive_trees.jpg width=277 height=200>
    <img src=Images/Symptoms_on_oleander.jpg    width=305 height=200>
</div>

The effects of the disease were first discovered and reported by Newton B. Pierce in the late 19th century.
The identification of the bacterium occurred only much later, in the 1970s.
Since then, the bacterium and the disease have been widely studied. 

The bacterium is typically found in areas ranging from tropical to temperate.
It is native to the Americas, but has also been detected in Europe, the Middle East and in the Far East. 
In the Americas, it is responsible for Pierce's disease in grapevines, citrus variegated chlorosis and coffee leaf scorch, while in Europe it is mainly responsible for **Olive Quick Decline Syndrome** (OQDS). 
In 2013 it was detected in the southern region of Apulia, in Italy, and has since spread all over the Salentine peninsula (the "heel" of the "boot").

There are multiple strains of the bacteria: _Xylella fastidiosa_ subsp. _fastidiosa_, subsp. _pauca_, subsp. _multiplex_, subsp. _sandyi_, etc.
The strain that is prevalent in Apulia is the **_pauca_** subspecies. 

The bacterium is spread by spittlebugs, froghoppers and sharpshooter leafhoppers, which are all xylem sap-feeding insects. 
**Philaenus spumarius**, the meadow froghopper or spittlebug, is considered the most important vector in Italy and in the EU, since it is the only one (proven so far) to transmit the bacterium in natural conditions.
These insects thrive and lay eggs in weeds in the spring months.
The young spittlebugs — nymphs — are born in the protected environment of a foam nest and, as summer approaches, they grow out of the dried out foam nest into adults.
The nymphs have reduced mobility compared to adult vectors (which can jump up to 70 cm), and feed only on herbaceous hosts (weeds), while the latter feed also on woody hosts.
Therefore, monitoring campaigns focus mostly, if not completely, on adult vectors. 
The figure below ([[2]](#2)) summarizes the vectors' life cycle. 

<div align="center">
  <img src=Images/Vector_lifecycle.jpg/ width=700 height=402>
</div>

There are hundreds of plants (belonging to different families) that are susceptible to the bacterium, according to data collected by the European Food Safety Authority (EFSA).
Not all plants are susceptible to all strains of the bacteria, and incubation periods/latency times (delay between infection time and onset of symptoms) vary across host species.
According to the [Food and Agriculture Organization](https://www.fao.org/in-action/saving-mediterranean-olives/en/#:~:text=Xylella%20fastidiosa%3A%20a%20spreading%20threat&text=The%20disease%20is%20difficult%20to,to%20more%20than%20a%20year.) (FAO), the incubation period of the bacteria can last from 7 months to more than a year, further complicating their detection and the implementation of preventive measures.



## Emergency measures

Since the initial detection of the _pauca_ strain of the bacterium near Gallipoli, in southern Salento, in 2013, the disease has caused the death of an estimated 21 million olive trees, from a total population of about 60 million.
Emergency measures implemented in Apulia vary according to the demarcated zones.
In the **infected zone** (_zona infetta_ in Italian), which includes the provinces of Lecce, Brindisi and part of Taranto, the disease cannot be eliminated and there is no mandate to eradicate infected plants. 
The **containment zone** (_zona di contenimento_) corresponds to the last 20 km of the infected zone.
Here, trees are constantly monitored: infected trees and all of its neighbors in a radius of 100 meters (including healthy trees) are eradicated.
Further measures include actively monitoring vector populations and maintaining a high level of surveillance over movement of plant materials and soil.
The **buffer zone** (_zona cuscinetto_) is an uninfected area that surrounds the containment zone for 10 km. 
In this area, authorities implement the same eradication and vector control measures of the containment zone, but they also take preventive steps such as removal of weeds, ploughing of the soil, diversified farming and planting of Xylella-resistant varieties of affected trees. 

<div align="center">
    <img src=Images/Demarcated_zones.jpg width=424 height=300>
</div>

In the first years of the outbreak, despite the EU drafting guidlines to combat the disease, Italian authorities were slow and inefficient in implementing them.
Most significantly, in 2015, charges were controversially filed against 10 individuals, including the Special Commissioner for the emergency and a number of head researchers at research institutions in the region, accusing them of aiding the propagation of the disease, through negligence and mismanagement.
The charges were dropped in 2019, when prosecutors failed to prove the causal link between the spreading of the disease and the actions of the investigated individuals. 



## Proposed model

## Data and analysis

## References
<a id="1">[1]</a>
White SM, Bullock JM, Hooftman DAP, Chapman DS. _Modelling the spread and control of Xylella fastidiosa in the early stages of invasion in Apulia, Italy._ Biol Invasions. 2017;19(6):1825-1837. doi: 10.1007/s10530-017-1393-5. Epub 2017 Feb 21. PMID: 32025190; PMCID: PMC6979717.

<a id="2">[2]</a>
EFSA (European Food Safety Authority), Vos S, Camilleri M, Diakaki M, Lázaro E, Parnell S, Schenk M, Schrader G and Vicent A, 2019. _Pest survey card on Xylella fastidiosa._ EFSA supporting publication 2019: 16(6):EN-1667. 53 pp. doi: 10.2903/sp.efsa.2019.EN-1667

<a id="3">[3]</a>
_Olive Oil Times_ (https://www.oliveoiltimes.com/)

<a id="4">[4]</a>
_European Commission_ (https://food.ec.europa.eu/plants/plant-health-and-biosecurity/legislation/control-measures/xylella-fastidiosa_en)

<a id="5">[5]</a>
_Emergenza Xylella (Apulia Region)_ (http://www.emergenzaxylella.it/)
