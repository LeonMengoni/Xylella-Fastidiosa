# Modeling the spread of Xylella Fastidiosa in Apulia, Italy

In this project we adopt the kernel spread model proposed by White et al. (2017, [[1]](#1)) in their study "Modelling the spread and control of Xylella fastidiosa in the early stages of invasion in Apulia, Italy" to estimate the speed of disease spread and to analyze the effect of the implemented emergency measures. 
Through this analysis, we aim to determine whether more effective strategies could have been employed in managing the outbreak.



## What is Xylella?

**_Xylella fastidiosa_** is a bacterial pathogen that is transmitted by insect vectors feeding on xylem sap from trees. 
Xylem cells act as tubes that transport water and nutrients (sap) from the soil up to the branches and leaves. 
The bacterium acts by obstructing the flow of nutrients through the xylem network, causing rapid dieback, starting from the leaves and proceeding up to the main trunk, resulting in the eventual death of the tree.
Common symptoms are leaf scorching, stunting and wilting, and a reduction in the size and quality of fruit production.
The figures below (taken from [[2]](#2)) capture the symptoms that detail the progress of the disease in the leaves (A), branches (B) and tree canopy (C) in the olive tree and oleander, respectively.  

<div align="center">
    <img src=Images/Symptoms_on_olive_trees.jpg width=277 height=200>
    <img src=Images/Symptoms_on_oleander.jpg    width=305 height=200>
</div>

The effects of the disease were first discovered and reported by Newton B. Pierce in the late 19th century.
The identification of the bacterium occurred only much later, in the 1970s.
Since then, the bacterium and the disease have been widely studied, albeit **no cure has been found**.

The bacterium is typically found in areas ranging from tropical to temperate.
It is native to the Americas, but has also been detected in Europe, the Middle East and in the Far East. 
In the Americas, it is responsible for Pierce's disease in grapevines, citrus variegated chlorosis and coffee leaf scorch, while in Europe it is mainly responsible for **Olive Quick Decline Syndrome** (OQDS). 
In 2013 it was detected in the southern region of Apulia, in Italy, and has since spread all over the Salento peninsula (the "heel" of the "boot").

There are multiple strains of the bacteria: _Xylella fastidiosa_ subsp. _fastidiosa_, subsp. _pauca_, subsp. _multiplex_, etc.
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



## Emergency measures in Apulia

The _pauca_ strain of the bacterium was first detected in 2013 in the province of Lecce, in Salento, but studies suggest that it was introduced between 2008 and 2010 through the port of Gallipoli (Kottelenberg et al, 2021).
Since its initial detection, the bacterium has caused the death of an estimated 21 million olive trees, out of a total population of about 60 million, resulting in a 50% reduction of olive oil production and an estimated economic loss of €1.6 billion ([[3]](#3)).
Emergency measures implemented in Apulia vary according to the demarcated zones.

In the **infected zone** (_zona infetta_ in Italian), which includes the provinces of Lecce, Brindisi and part of Taranto, the disease cannot be eliminated and there is no mandate to eradicate infected plants. 

The **containment zone** (_zona di contenimento_) corresponds to the last 20 km of the infected zone.
Here, trees are constantly monitored: infected trees and all of its neighbors in a radius of 100 meters (including healthy trees) are eradicated.
Further measures include actively monitoring vector populations and maintaining a high level of surveillance over movement of plant materials and soil.

The **buffer zone** (_zona cuscinetto_) is an uninfected area that surrounds the containment zone for 10 km. 
In this area, authorities implement the same eradication and vector control measures of the containment zone, but they also take preventive steps such as removal of weeds, ploughing of the soil, diversified farming and planting of Xylella-resistant varieties of affected trees. ([[4]](#4), ([[5]](#5))

<div align="center">
    <img src=Images/Demarcated_zones.jpg width=424 height=300>
</div>

In the first years of the outbreak, despite the EU guidelines, Italian authorities were slow and inefficient in implementing control measures. 
Most significantly, in 2015, following widespread protests led by local farmers and enviornmental organizations against the felling of monumental olive trees, eradication measures were halted.
Controversially, charges were filed against 10 individuals, including the Special Commissioner for the emergency Giuseppe Silletti and Donato Boscia, one of the discoverers of Xylella in Italy and lead researcher of the National Research Council (CNR), accusing them of aiding the propagation of the disease, through negligence and mismanagement.
The charges were dropped in 2019, when prosecutors failed to prove the causal link between the spreading of the disease and the actions of the investigated individuals, but by then, the damage was done and the disease had spread almost uncontrolled in Salento, the lower half of the region. 

Recently, scientific authorities have finally acknowledged that the implemented measures are bearing fruit.
Not only has the northward spread of the disease slowed down, but disease symptoms seem to be milder, and the incidence of Xylella-carrying vectors has also decreased.
Eradication of Xylella is impossible, therefore the only solution is **coesixtence**. 
With this in mind, research is focusing on:

- developing more efficient ways of detecting the disease, by using dogs, drones or DNA tests
- finding and planting Xylella-resistant (if not immune) varieties of tree species
- developing innovative and more sustainable vector control practices.

## Proposed model

To study the spatial spread of the disease, we adopted the model proposed by White et al. in their 2017 paper ([[1]](#1)).

The model runs over a 1km² resolution grid of the Apulia region, with a temporal scale of 1 year, which accounts for the seasonality of the vector's lifecycle. 
Spreading occurs in two phases: local growth within a cell and spatial spread between cells. 
Control mesasures are implemented only in the control zone. 

In practice, in the code we apply the different steps consecutively, therefore updating the infectious population step-by-step.
In other words, for every timestep $t$: 

$$
I_t \rightarrow I^G_{t+1} \rightarrow I^S_{t+1} \rightarrow I^L_{t+1} \rightarrow
\begin{cases}
I_{t+1} & \text{if no control measures}\\
I^C_{t+1} \rightarrow I_{t+1} & \text{if control measures}
\end{cases}
$$

where $G$ recalls the local growth, $S$ the short-distance dispersal, $L$ the long-distance dispersal and $C$ the application of control measures. 

### 1. Local growth

Local growth of the density of infected trees is modeled by the Gompertz function $I^G(t) = K \textrm{e}^{-B\textrm{e}^{-At}}$, which, when discretized on an annual time scale, is written as:

$$
I^G_{t+1}(x,y) = K(x,y) \left(\frac{I_t(x,y)}{K(x,y)}\right)^{\textrm{e}^{-A}} = f(I_t(x,y))
$$

where $A$ is the rate of population growth, $B$ is related to the initial proportion of plants that are infected, and $K(x,y) = d(x,y) + a(1 - d(x,y))$ is the carrying capacity; $d(x,y)$ is the proportional cover/density of olive groves in a 1km² grid cell, and $a\in[0,1]$ is the carrying capacity in non-olive grove habitat, relative to that in olive groves. 
Here, $I_t(x,y)$ is the density of infected trees, such that $I_t(x,y)\in[0, d(x,y)]$. 
The incidence is therefore defined as: 

$$
i_t(x,y) = \frac{I_t(x,y)}{d(x,y)}
$$

The Gompertz function parameters have been fitted on the data, yielding the values $\lbrace A = 3, B = 14.069\rbrace$, while we assume $a = 0$ (we consider only olive trees, no other trees) and $d(x,y)$ is taken from the file ```olivegrowthprop.mat```.

### 2. Spatial spread

Spatial spread of the disease follows a _stratified dispersal_ pattern, i.e. a two-process dispersal involving a short-distance dispersal, which represents the local mobility of jumping vectors, and a long-distance dispersal, which accounts for unintentional flight dispersal of vectors by wind or hitchhiking in cars and trucks.

#### Short-distance dispersal
The short-distance dispersal is modeled by a deterministic 2D kernel:

$$
I^S_{t+1}(x,y) = \sum_i\sum_j\hat{k}(x-i,y-j)I^G_{t+1}(i,j) \qquad \Longrightarrow \qquad I^S_{t+1} = \hat{k} * I^G_{t+1}
$$

where the sum is clearly over the dimensions of the grid. 

We have tried different kernels, namely an exponential and a gaussian:

$$
\begin{aligned}
& \hat{k}_e(x,y) = \textrm{e}^{-\frac{(x^2 + y^2)^{\frac{1}{2}}}{\beta}} \\
& \hat{k}_g(x,y) = \textrm{e}^{-\frac{x^2 + y^2}{2\beta^2}}
\end{aligned}
$$

where the mean dispersal distance $\beta$ (measured in kms) is assumed to be $\lbrace\beta = 0.1\rbrace$.
The kernel is not normalized because we want the population in the origin cell $(x,y)$ to remain the same. 
If we were modeling the movement of infected individuals that spread across the map, then we would have to normalize the kernel.

#### Long-distance dispersal

The long-distance dispersal is modelled by isotropic stochastic dispersal. 
In other words, every cell has a probability given by $u(x,y) = \rho(x,y) I_t(x,y)$ (where $\rho(x,y)\sim\mathcal{U}[0,1]$) to generate a random number of dispersers. 
If $u(x,y) > p$ ($p$ sets a threshold probability), then the cell $(x,y)$ randomly generates $M\leq M_{max}$ dispersers, which disperse a random distance drawn from a 2D gaussian distribution $\mathcal{N}(0,D)$ ($[D] = \textrm{km}$).
Finally, a cell that is reached by one of the random dispersers is further infected according to the initial proportion of infected according to the Gompertz function:

$$
\Delta I_t(x,y) = (d(x,y) - I_t(x,y))\textrm{e}^{-B}
$$

The parameters are set as $\lbrace p = 0.2, M_{max} = 5, D = 20\rbrace$.

If we wanted to model the average infection level ($\overline{I}_t$) over $N\rightarrow\infty$ runs of the simulation, we can find an analytical expression for the long-distance kernel. 
The evolution is written as the following (I keep $I_t$ instead of $\overline{I}_t$ just because I'm too lazy to go change all the equations): 

$$
\begin{aligned}
I_{t+1}(x,y) & = I_t(x,y) + \Delta I_t(x,y) \\
& = I_t(x,y) + \overline{M}_{in}(x,y)\cdot(d(x,y) - I_t(x,y))e^{-B}
\end{aligned}
$$

where $\overline{M}_{in}(x,y)$ is the average number of dispersers jumping into the cell $(x,y)$. 
This term can be written in the following way:

$$
\overline{M}\_{in}(x,y) = \sum\_i\sum\_j p(x-i,y-j) \overline{M}_{out}(i,j)
$$

where the probability $p(x-i,y-j)$ follows a discretized gaussian distribution, .
The final term $\overline{M}_{out}(i,j)$ is the average number of dispersers jumping out of cell $(i,j)$ and is equal to the average number of dispersers of $(i,j)$ (if it disperses) multiplied by the probability that $(i,j)$ disperses: 

$$
\begin{aligned}
\overline{M}\_{out}(x,y) & = \frac{M_{max}}{2}\cdot P(\rho I_t(x,y) > p) \\
& = \frac{M_{max}}{2}\cdot (1 - P(\rho I_t(x,y) < p)) \\
& = \frac{M_{max}}{2} \left[1 - \frac{p}{d(x,y)}\left(\ - \text{log}\frac{p}{d(x,y)}\right)\right]
\end{aligned}
$$

The expression 

$$
P(\rho I_t(x,y) < p) = 1 - \frac{p}{d(x,y)}\left(1 - \text{log}\frac{p}{d(x,y)}\right)
$$ 

is the cumulative distribution of the product of two uniformly distributed random variables, respectively $\rho$ on $\[0,1\]$ and $I_t(x,y)$ on $\[0,d(x,y)\]$ for all times $t$.
The latter is a quite strong assumption considering that, as the infection progresses, the fraction of cells with a high level of infection ($`I_t(x,y)\lesssim d(x,y)`$) increase (so the distribution of $I\_t(x,y)$ would be skewed towards $`d(x,y)`$).
However, in areas where the infection has progressed, random dispersers would mainly disperse into cells with a high level of infection, and their contribution would be negligible to the overall spreading of the infection.
Therefore, we practically imagine that the uniform distribution of $I_t(x,y)$ is applied only to the front of the infection.

To summarize, the final expression for the average infections $\overline{I}^L\_{t+1}=I^L_{t+1}$ is:

$$
\begin{aligned}
I^L_{t+1}(x,y) & = I^S_{t+1}(x,y) + \Delta I^S_{t+1}(x,y) \\
& = I^S_{t+1}(x,y) + \overline{M}\_{in}(x,y)\cdot(d(x,y) - I^S\_{t+1}(x,y))e^{-B} \\
& = I^S_{t+1}(x,y) + (d(x,y) - I^S\_{t+1}(x,y))e^{-B}\cdot\sum\_i\sum\_j p(x-i,y-j) \overline{M}\_{out}(i,j) \\
& = I^S_{t+1}(x,y) + (d(x,y) - I^S\_{t+1}(x,y))e^{-B}\cdot\sum\_i\sum\_j p(x-i,y-j) \frac{M_{max}}{2} \left[1 - \frac{p}{d(i,j)}\left(1 - \text{log}\frac{p}{d(i,j)}\right)\right] \\
& = I^S_{t+1}(x,y) + (d(x,y) - I^S\_{t+1}(x,y))e^{-B}\frac{M_{max}}{2}(p(x,y) * g(d(x,y)))
\end{aligned}
$$

where $g(d(x,y)) = \left[1 - \frac{p}{d(x,y)}\left(\ - \text{log}\frac{p}{d(x,y)}\right)\right]$.

In absence of control measures, $I_{t+1} = I^L_{t+1}$.

### 3. Control measures

In the model, the demarcated areas are subdivided in a slightly different way. We have four zones: infected (**IZ**), eradication (**EZ**), buffer (**BZ**) and surveillance (**SZ**). 
Control measures are implemented only in the eradication and buffer zones, which together form the control zone (**CZ**). 
In this zone, we assign a probability of infection detection $p_{detect}(x,y)\sim\mathcal{U}[0,1]$ to every cell.
If the surveillance efficiency $s\in[0,1]$ (given as a parameter) is greater than the detection probability ($s > p_{detect}$), the infected trees in the cell are eradicated and not replaced.

$$
I^C_{t+1}(x,y) = 
\begin{cases}
0 & \text{if} p_{detect}(x,y) < s \\
I^L_{t+1}(x,y) & \text{otherwise}
\end{cases}
$$

On average, we would obtain 

$$
\overline{I}^C\_{t+1} = (1 - s)\overline{I}^L_{t+1}
$$

The surveillance efficiency in the two zones within the control zone may vary: for example, it can be higher in the eradication zone and lower in the buffer zone.   

The parameter here is therefore only $\lbrace s\rbrace$.


## Data and analysis

### 1. Estimation of bacteria spread speed
In a paper by Kottelebnerg, Saponari (2021), the speed was estimated as 10.0 km per year (95% confidence interval: 7.5–12.5 km per year). 

We use as method a simple linear regression over maximum distances from the epicentre over time. 
We get an estimate for the speed over all years (2013-2023) or only over the first 4 years of samples (2013-2016).

An alternative method could be to use Kottelennerg's method with a logistic function. 
### 2. Obtain parameters from least squares method by comparing data speed to simulation speed
### 3. Introduction of latency period 
In White (2020), the incubation period (infected but asymptomatic, and negligible to no infectivity) is estimated with Bayesian methods at 1.2 years (1-1.3 95% credibility interval).

**Definition of risk**

## Summary of other main results and findings 

White et al, 2017
Kottelenberg et al, 2021

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
