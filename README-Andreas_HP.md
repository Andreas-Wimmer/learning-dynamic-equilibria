# Learning-dynamic-equilibria
 This repository contains all the code that I need for my master thesis.

 ## The basic idea
 The basic idea is that we have a dynamic s-t-netowrk with travel times and capacities given, where s is some origin/source node/point and t is some destination/sink node/point. At the source/origin players, which we will call agents, start and they choose a s-t-path along which they traverse the network and they experience costs in form of travel time, partially because of the travel times of the edges and partially because of queues that build up along the path they chose. In a time interval $[0, T]$ that could for example mimic a whole day, there is now a whole population of agents $[0,r]$ that starts at the origin and have to choose a path. Once this is done, on "the next day", every agent reevaluates his choice and maybe takes a different path this "day", because it was shorter "yesterday". The dynamics that evolves out of that process is a learning dynamics. The overall question now is, if this learning dynamics converges and if yes, if the resulting state is stable in the sense that no agent has an incentive to deviate. Such a stable state is then called dynamic equilibrium.

 ## The rough algorithm

 To implement that, one has to do the following steps:

 1. Initialize the flow in some way, for example according to the minimal capacitiy on each path.
 2. Compute the regularized best response to the so called average population frequency defined as $\bar{X}_n = \frac{1}{n-1} \cdot \sum_{i=1}^{n-1} X_i$, where $X_i$ refers to the decision profile on the $i$-th "day". The regularized best-response is in the discretization defined as the solution of the following minimization problem:
 $ \min\{\sum_{p \in \mathcal{P}} \sum_{j \in J} \Psi_p(\bar{X}_n)[j]\cdot h_p[j] + \epsilon \cdot (h_p[j] - (\bar{X}_n)_p[j])^2 \vert h_p[j] \geq 0 \forall p \in \mathcal{P}, j \in J, \sum_{p \in \mathcal{P}} h_p[j] = r \forall j \in J\}$,
 where $J$ is a set of intervals where we assume the wante flow is constant on and $\Psi$ is a path delay operator and $\varepsilon$ is a regularization weight.
 3. Update the current flow and the population average and compare the current flow with the last flow.
 4. Repeat steps 2&3 until either one has reached a given number of iterations or the comparison in step 3 yields that one has reached a given accuracy of the convergence.

## Special cases
There are two special cases, namely when you either have only one congested edge per path or only one path per congested edge, where we know theoretically that the above procedure gets us to a dynamic equilibrium or in practice at least to an $\varepsilon$-regularized dynamic equilibrium.

## General networks
We know that in general networks we do not get the same properties as in the above described special cases and hence also do not know how the above described learning dynamics behaves. 
