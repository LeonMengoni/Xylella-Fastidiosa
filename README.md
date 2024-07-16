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
    <img src=Images/Yearly_samples_2013-2023.png width=1000 height=844>
</div>

#### i. Linear regression

In this first method, we use a simple linear regression over the evolution of the maximum distance from the epicentre over time. 
We get an estimate for the speed over all years (2013-2023) or only over the first 4 years of samples (2013-2016).
This method is, in truth, very rough and subject to bias.


The analysis of positive samples yielded the following results: 

$$
\begin{aligned}
v & = 9.5\pm0.4 \text{km}/\text{year} \\
v_{\text{init}} & = 16.3\pm0.5 \text{km}/\text{year}
\end{aligned}
$$

<div align="center">
    <img src=Images/Disease_Spreading_Speed.png width=700 height=417>
</div>

#### ii. Kottelenberg method

In a paper by Kottelebnerg, Saponari (2021), the speed was estimated as 10.0 km per year (95% confidence interval: 7.5–12.5 km per year).

Apply a logistic function to the shape of the epidemic front, assuming that the disease spread speed remains constant over time 

### 2. Obtain parameters from least squares method by comparing data speed to simulation speed

#### Sensitivity analysis

### 3. Introduction of latency period 
In White (2020), the incubation period (infected but asymptomatic, and negligible to no infectivity) is estimated with Bayesian methods at 1.2 years (1-1.3 95% credibility interval).

**Definition of risk**

## Summary of other main results and findings 

White et al, 2017
Kottelenberg et al, 2021

## Extra
### Long-distance dispersal equation
If we wanted to model the average infection level ($\overline{I}_t$) over $N\rightarrow\infty$ runs of the simulation, we can find an analytical expression for the long-distance kernel. 
The evolution is written as the following: 

$$
\begin{aligned}
\overline{I}\_{t+1}(x,y) & = \overline{I}\_t(x,y) + \Delta \overline{I}\_t(x,y) \\
& = \overline{I}\_t(x,y) + \overline{M}_{in}(x,y)\cdot(d(x,y) - \overline{I}_t(x,y))e^{-B}
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
& = \frac{M_{max}}{2} \left[1 - \frac{p}{d(x,y)}\left(1 - \text{log}\frac{p}{d(x,y)}\right)\right]
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
