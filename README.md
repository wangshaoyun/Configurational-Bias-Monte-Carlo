
![Schmetical Diagram](https://github.com/wangshaoyun/Configurational-Bias-Monte-Carlo/blob/master/cylinder.png "Simulation System")
# Adsorption of Weak Polyelectrolytes
## About the program
1. This is a program used to simulate adsorption of weak polyelectrolytes (PE) on cyliner, plate or sphere.
2. Configuration biased Monte Carlo method [1] is used to acceralate equillibrium. There are 4 moving forms to improve the sampling efficience in middle monomers [2]. The first is retrace part chains from start point and add them to start point. The second is retrace part chains from start point and add them to the end point. The third is retrace part chains from end point and add them to end point. The fourth is retrace part chains from end point and add them to the start point.
3. The interactions include Lennard-Jones potential, finite extensile nonlinear elastic (FENE) potential, bending stiffiness, torsion stiffiness and Coulomb potential [3].
4. Constant pH titration method is used to simulate the ionization process.
5. The short-range potential is calculated by pure cell lists method [4].
6. Multiple time step method [5] which separate the Coulomb potential into short-range part and long-range part. Long-range part is updated after several Monte Carlo steps and is calculated by smooth particle mesh Ewald method [6].
>[1] Frenkel D, Klein M, Parrrinello M. "Understanding Molecular Simulation: From Algorithm to Applications". 2nd ed. Singapore: Elsevier Pte Ltd, 2002.  
>[2] Allen M P, Tildesley D J. Computer simulation of liquids. 2nd ed. Oxford: Oxford University Press, 2017.  
>[3] Zhang F, Wang S, Ding H, Tong C. Simulations of 3-arm polyelectrolyte star brushes under external electric fields. *Soft Matter*, **15** (12), 2560-2570, 2019.  
>[4] S. Y. Wang, C. H. Tong. Cell-lists Method for Monte Carlo Simulation, to be submitted.  
>[5] K. Bernacki, B. Hetenyi, B. J. Berne. "Multiple "Time Step" Monte Carlo Simulations: Application to charged systems with Ewald Summation." *Journal of Chemical Physics*, **121** (1), pp.44-50, 2004.  
>[6] Ulrich Essmann, Lalith Perera, Max L. Berkowitz, Tom Darden, Hsing Lee, Lee G Pedersen. "A Smooth Particle Mesh Ewald Method." *Journal of Chemical Physics*, **103** (19), pp. 8577-8593, 1995.  

## About Parameter Files 
+ energy_data.txt: It is a parameter file which includes Bjerrum length, external electric field
+ system_data.txt: It is a parameter file which includes system, time, brushes parameter.  

## Compiling files
1. Determine the parameters such as box size, brush chians, time step and so on.
2. Set the parameters in energy_data.txt and system_data.txt
3. Open terminal and enter the current file content
4. To compile the files:
```
$ make
```
  
5. To execute the files:
```
$ ./main &
```


