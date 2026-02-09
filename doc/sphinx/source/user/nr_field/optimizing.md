# Optimising Nr Manually
To rigorously optimise the Nr field of a simulation, [wavefield learning](learning.md) has to be performed during a full resolution run. But, in practice, users develop a feeling for the sort of Nr values to expect with the types of models they simulate. They are then able to save computational cost from the get go by capping Nr to an expected maximum value. Unfortunately, at this stage, this is really just a matter of experience and we cannot provide any reliable rules of thumb. But we can supply some practical tips from our own experience:

* **Basic:** the standard approach is to use a Nr bound by inplane resolution and a constant cap, i.e.:
    * *type_Nr: CONSTANT*
    * *constant: an educated guess for the maximum possible Nr-value*
    * *bound_Nr_by_inplane: true*
* **Advanced:** for global models, depth-dependent variations of Nr can be predicted rather well (shown in Leng et al. 2019), so it may be better to use the analytical field rather than a constant cap Nr, i.e.:
    * *type_Nr: ANALYTICAL*
    * *control_depths: vector of depths where max Nr is known (unit in metres, vector given as [depth1, depth2, ...])*
    * *Nr_at_control_depths: educated guess for the maximum possible Nr-value at the above depths (also given as vector)*
    * *bound_Nr_by_inplane: true*
* Note: the analytical section is specifically written to allow for easy additional user input in case they observe other predictable variations
* **Ballpark:** Hard to say as this depends on model smoothness, size, the source frequency and many other factors. But to give some numbers, below 500 is considered a rather small run and above 1500 is considered rather big.
* **Always double check:** It comes at little extra cost to switch on wavefield learning when performing a simulation. This allows you to double check your initial guess for the Nr field and, if everything worked well, you may even use the resulting map to fully optimise your run in the future. Details on what to look our for when checking your initial Nr field are given [here](learning.md).
* **Changing source frequency:** We have fared well with the assumption that doubling the source frequency also doubles the maximum Nr
