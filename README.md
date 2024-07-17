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
The vectors can sometimes also spread long-distances, by hitchhiking on cars or trucks. 
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

Depending on the value of $A$, the incidence of a cell grows more or less rapidly, even for very small initial seeding.

<div align="center">
  <img src=Images/Local_growth_A.png width=400 height=400>
</div>

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
In absence of control measures, $I_{t+1} = I^L_{t+1}$.

### 3. Control measures

In the model, the demarcated areas are subdivided in a slightly different way. We have four zones: infected (**IZ**), eradication (**EZ**), buffer (**BZ**) and surveillance (**SZ**). 
Control measures are implemented only in the eradication and buffer zones, which together form the control zone (**CZ**). 
In this zone, we assign a probability of infection detection $p_{detect}(x,y)\sim\mathcal{U}[0,1]$ to every cell.
If the surveillance efficiency $s\in[0,1]$ (given as a parameter) is greater than the detection probability ($s > p_{detect}$), the infected trees in the cell are eradicated and not replaced.

$$
I^C_{t+1}(x,y) = 
\begin{cases}
0 & \text{if } p_{detect}(x,y) < s \\
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

Firstly, we attempt to calculate the spreading speed of the disease from the monitoring data. 
On the model, we perform a sensitivity analysis to identify the key parameters that cause this speed to change most significantly.
Then, we search for the model parameters that most closely resemble our speed estimates.
**Finally, we introduce control measures to explore variations in disease spreading speed.**

### 1. Estimation of disease spreading speed from monitoring data

We attempt two methods to estimate the disease spreading speed: a linear regression over maximum distances from the epicentre and a logistic fit to the epidemic front. Both methods are limited in the attendability of their results.

The initial moment of introduction of the bacterium is unknown, as well as its epicentre.
The first positive samples were taken in the _comune_ of Trepuzzi, in the hinterland of the province of Lecce, but pretty soon samples were also taken in and around Gallipoli.
Subsequent studies (Kottelenberg, 2021) have shown that Gallipoli is the most likely point of introduction of the bacterium, in the timeframe 2008-2010.
In addition, monitoring data is extremely heterogeneous, both in time and in space.
As can be seen in the figure below, sampling areas change drastically over the years: while sampling occurred more uniformly over the whole region in the first couple of years, after 2015 most samples have been taken only in the buffer and containment zones, and close to none in the infected zone.
Also, the number of samples taken every year has fluctuated a lot.

In the following, the epicentre of the disease is assumed to be Gallipoli (40.055851°N, 17.992615°E).

<div align="center">
    <img src=Images/All_samples_by_year.png width=1000 height=750>
</div>

#### _i_) Linear regression

In this first method, we use a simple linear regression over the evolution of the maximum distance from the epicentre over time. 
We get an estimate for the speed $c$ over all years (2013-2023) or only over the first 4 years of samples (2013-2016), $c_{\text{init}}$.
This method is, in truth, very rough and subject to bias, for the reasons stated above, regarding the missing data. 

The analysis of positive samples yielded the following spreading speeds: 

$$
\begin{aligned}
c & = 9.5\pm0.4 \text{km}/\text{year} \\
c_{\text{init}} & = 16.3\pm0.5 \text{km}/\text{year}
\end{aligned}
$$

<div align="center">
    <img src=Images/Disease_Spreading_Speed.png width=700 height=417>
</div>

#### _ii_) Kottelenberg method

For the second method, we take inspiration from Kottelenberg et al. (2021, [[6]](#6)).
Here, the authors attempt to estimate the shape of the invasion front and its rate of movement. 

The monitoring data was grouped in 1-km wide distance classes from the epicentre of Gallipoli and we calculated the proportion of positive tests in each class. 
We chose for the shape of the disease front a deterministic function, namely a logistic function of the type:

$$
p = p(x,t) = \frac{1}{1 + e^{x-(x_{50}) + ct}}
$$

where the function depends also on time $t$ (years elapsed from start of epidemic, i.e. 2013), and where $c$ is the spreading speed of the disease (and is assumed to be constant through the years) and $x_{50}$ is the (negative) x-value (distance from Gallipoli) of the half-maximum of the curve at $t = 0$. 

In each distance class $d$, the number $n_d$ of samples tested is known, while the number of positives $pos_d$ is a binary stochastic variable.
Therefore, we choose the binomial distribution for fitting the model to the data.
Finally, we evaluate the model parameters by maximizing the likelihood (minimizing the negative log-likelihood).

To summarize, we proceed as follows. 
We only analyze the first 6 years, 2013-2018, because data becomes very clustered in the buffer zone for later years.
For every one of these years, we minimize the negative log-likelihood to obtain the fit parameters: 

$$
\begin{aligned}
NLL & = -\sum_{d\in\lbrace 1,2,\cdots,d_{max}\rbrace} \log[\mathcal{B}(pos_d, n_d)] = -\sum_d \log\left[{n_d\choose pos_d} p^{pos_d} (1-p)^{n_d-pos_d} \right]
\end{aligned}
$$

<div align="center">
    <img src=Images/Speed_from_data_fit.png width=1000 height=667>
</div>

The parameter we are most interested in is $c$:

$$
\begin{aligned}
c_{2013} & = 0.2\text{km}/\text{year} \\
c_{2014} & = 8.7\text{km}/\text{year} \\
c_{2015} & = 6.1\text{km}/\text{year} \\
c_{2016} & = 5.9\text{km}/\text{year} \\
c_{2017} & = 6.4\text{km}/\text{year} \\
c_{2018} & = 5.2\text{km}/\text{year} 
\end{aligned}
$$

As we can see, the estimate for $c$ in 2013 is extremely off, due to the poor data.
For the following years, estimates for $c$ are in a reasonable range.

In the paper, a more sophisticated analysis led to an estimate of $c_{Kott.} = 10.0\text{km}/\text{year}$ (with 95% confidence interval $7.5-12.5\text{km}/\text{year}$).  

### 2. Sensitivity analysis on simulation model 

To perform a sensitivity analysis for the model, we study how modifying the input parameters affect the model outputs (i.e. the disease spreading speed $c$).
We therefore have to also calculate the speed $c_{\text{sim}}$ from our simulation.

#### Estimation of disease spreading speed from simulation

We propose to calculate the speed differently from the case of actual data, due to the fact that in the simulation we have all the "data".
Instead of using the maximum distance from the epicentre, we use the average distance of all infected cells, and perform a linear regression with respect to time.

$$
c_{\text{sim}} = \frac{1}{|\text{infected cells}|} \sum_{(x,y)\in\lbrace\text{infected cells}\rbrace} \frac{d_O(x,y; t)}{t}
$$

where $d_O(x,y;t)$ is the distance of cell $(x,y)$ from the epicentre $O$, at time $t$.

To eliminate (however possible) the stochasticity of the model, we define the **risk** $R(x,y;t)$ (for every cell $(x,y)$ and every timestep $t$) as the average incidence over $N$ runs of the simulation, where $N$ is typically in the order of ~ 10.
(White et al. [[1]](#1) define it over $N=10000$ runs, but we will keep this number low because of constraints on time and computational resources.)
When calculating the risk, compared to a single run of the simulation, there will be cells that are further out with a positive average incidence, therefore skewing the estimate of the speed $c_{\text{sim}}$ towards higher values.
To mitigate this effect, we weight the distance of a cell from the epicentre by its average incidence (risk), as follows:

$$
c_{\text{risk}} = \frac{1}{|\text{infected cells}|} \sum_{(x,y)\in\lbrace\text{infected cells}\rbrace} \frac{R(x,y;t)\cdot d_O(x,y;t)}{t}
$$

#### Sensitivity analysis

We study how the disease spreading speed $c_{\text{risk}}$ varies as we tweak the model input parameters.
The parameters that we will explore are: $\lbrace A, B, \beta, p, M_{max}, D\rbrace$.
For evaluating the risk, we fix $N = 30$. 

##### _i_) $A$: rate of local growth

From the below figure we can see that increasing the rate of local growth $A$ has an effect on increasing the speed, up to a point, after which the speed doesn't vary much anymore. 
If we recall the graph showing the local growth as a function of $A$, we can easily understand why: for values of $A$ exceeding 4 or 5, one timestep of the simulation (i.e. one year) is enough to send the incidence close to 100%, and for larger $A$ there is not much change. 

<div align="center">
    <img src=Images/Speed_A.png width=600 height=484>
</div>

##### _ii_) $B$: (related to) initial proportion of infected

Let's recall that the parameter $B$ is related to the seeding of a new cell through the proportion $e^{-B}$.
Therefore, for greater $B$ we would obtain a smaller seeding.
Since we use a tolerance threshold to eliminate numerical noise, set at a low value of $1e^{-8}$, for very large values of $B$, new cells wouldn't get infected at successive iterations of the simulation.
This explains the effect of the dipping curve for large $B$. 
On the other hand, we see that for smaller values of $B$, the speed doesn't change all that much.
This is obvious since we just concluded in the previous section that it is $A$ that mainly drives the local growth.

<div align="center">
    <img src=Images/Speed_B.png width=600 height=450>
</div>

##### _ii(bis)_) $A$ and $B$ (simultaneously)

To further prove the previous statement, we plot the speed as a function of both parameters $A$ and $B$ and we see that the greater changes occur along the $A$-axis, and not the $B$-axis.

<div align="center">
    <img src=Images/Speed_AB.png width=600 height=600>
</div>

##### _iii_) $\beta$: mean dispersal distance

As we could imagine, by increasing the mean dispersal distance of the short-distance kernel, the speed of the epidemic increases, approximately linearly. 
Looking at the y-axis, we also see that this parameter is the greatest driver of the disease spreading speed. 

**We could try seeing the balancing between short-distance and long-distance dispersal distances, $\beta$ and $D$.**

<div align="center">
    <img src=Images/Speed_beta.png width=600 height=454>
</div>

##### _iv_) $M_{max}$: maximum number of long-distance dispersers per cell

Increasing the maximum number of dispersers $M_{max}$, the average speed increases, as can be expected. 

<div align="center">
    <img src=Images/Speed_M_max_30.png width=600 height=450>
</div>


##### _v_) $D$: variance of long-distance dispersal

Analogously, for $D$, the variance of the gaussian long-distance jumps of dispersers.

<div align="center">
    <img src=Images/Speed_D.png width=600 height=450>
</div>

##### _v(bis)_) $M\_{max}$ and $D$ (simultaneously)

Again, here we notice how the greater change is given by the variation of $D$, although, admittedly, we would have to explore the space of the parameters $M_{max}$ and $D$ a bit more to give a more definite answer. 

<div align="center">
    <img src=Images/Speed_M_maxD.png width=600 height=600>
</div>

### 3. Model parameter search

Now, we want to find a set of parameters that may give us back the values for the speed that we estimated for the data, by a least squares minimization.
We assumed a constant speed for the epidemic front, unchanged in every year; thus, we can use as estimate the average speed calculated from the years 2014-2018: $\overline{c} = 6.5\text{km}/\text{year}$.

Having called our model $\mathcal{M}(A, B, \beta, M_{max}, D)$ our objective is to find the set of parameters that minimizes:

$$
f(A, B, \beta, M_{max}, D | \overline{c}) = [c(\mathcal{M}(A, B, \beta, M_{max}, D)) - \overline{c}]^2
$$

To minimize the search time, we can eliminate $B$ from the analysis since we determined that it doesn't affect the final speed estimates very much.
We do not have an objective function to calculate $c(\mathcal{M}(A, B, \beta, M_{max}, D))$, so the gradient is not known.
Gradient methods for minimization cannot therefore be used, so we will use a 

Genetic algorithms could also be used for more efficient optimization. 

### 4. Varying control measures

Determine whether implemented measures have yielded results.

### 5. Further proposals: introduction of latency period 
In White (2020), the incubation period (infected but asymptomatic, and negligible to no infectivity) is estimated with Bayesian methods at 1.2 years (1-1.3 95% credibility interval).
One idea for further analysis would be to incorporate this information in the simulation model.
At present (17/07/2024), this was not implemented due to lack of time.

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

<a id="6">[6]</a>
Kottelenberg, D., Hemerik, L., Saponari, M. et al. Shape and rate of movement of the invasion front of Xylella fastidiosa spp. pauca in Puglia. Sci Rep 11, 1061 (2021).
