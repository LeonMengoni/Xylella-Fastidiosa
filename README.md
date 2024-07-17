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
    <img src=Images/All_samples_by_year.png width=1000 height=844>
</div>

#### i. Linear regression

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

#### ii. Kottelenberg method

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

### 2. Sensitivity analysis

### 3. Model parameter search 

**Definition of risk**

### 4. Further proposals: introduction of latency period 
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
