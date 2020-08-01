 
 
 
 ### Setup
 
 Please download and install Amber20. At the end, you should have a path to amber and set the environment variable:
 ```
export AMBERHOME=/Users/austin/Downloads/amber20
```
where inside the folder there is a file ```amber.sh```.


 
Create a new Conda env with python 3.7 and install the following packages below.
```
conda create -n myenv python=3.7 -y
conda activate myenv
conda install -c conda-forge -c omnia -c openeye openmm openeye-toolkits pymbar parmed ambertools
``` 


### Usage
There will be more features added in the future. For now, this is the code on the master branch. 

Please make sure your working directory is the base of this repo.

```
cd path-to-mmgbsa_scripts/
```

There are two phases to MMGBSA/PBSA calculations: equilibration (or relaxation as I think of it), and the production simulation
where we sample states to use for the MMGBSA calculation. MMGBSA calculations are computed on static frames, and averaed over those frames.
There are two approaches implemented. ```--calcFrames 5``` for example will sample five equally spaced frames from the production 
trajectories (this is the naive approach), ```--mbar 100``` for example will sample 100 equall spaced frames, and then use a technique
from John Chodera's lab to subsample those frames to "maximizes number of effectively uncorrelated samples" [1].  For instance,
```--mbar 100``` will sample 100 frames, and then subsample 20 frames that are the least correlated. By default, the code will use 
```--mbar``` instead of naive sampling. The default settings should be perfect for most cases. 


Here is an example. You can use -h to see more options if desired. 
PBSA is very slow. If you use ```--method pbsa``` please also do not use ```--mbar```, instead pass ```--calcFrames ``` with a number
no greater than 10 unless you have a lot of time on your hands... 

```
export AMBERHOME=...../amber20
export OE_LICENSE=...../oe_license.txt
source $AMBERHOME/amber.sh
conda activate myenv #do this after the source.
CUDA_VISIBLE_DEVICES=0 python run_mmgpbsa.py --pdb /homes/aclyde/MPro_0387_Gen3L_5.pdb -v 0 \ 
               --odir testdir --platform CUDA
```

If the script is being wonky, try passing ```--amber_path $AMBERHOME```, though it should automatically detect this from the env. 

#### Example output
```
/Users/austin/anaconda3/envs/mmpbsa_scripts/bin/python /Users/austin/PycharmProjects/mmpbsa_scripts/run_mmgpbsa.py --pdb /Users/austin/MPro-docked/MPro_0387_Gen3L_5.pdb --odir /Users/austin/PycharmProjects/mmpbsa_scripts/test --platform OpenCL
INFO [SystemLoader:split_complex_from_system] Reading input PDB /Users/austin/MPro-docked/MPro_0387_Gen3L_5.pdb
INFO [SystemLoader:split_complex_from_system] Split complex. atom sizes-- lig: 18, prot: 2370, water: 22, other: 331
INFO [SystemLoader:prepare_simulation] Built simulation using platform OpenCL with properties {'Precision': 'mixed'}
INFO [SystemLoader:prepare_simulation] Minimizing and setting velocities to 310.15 K
INFO [SystemLoader:prepare_simulation] Equilibrating for 5000 steps, or 10 ps
#"Progress (%)"	"Step"	"Potential Energy (kJ/mole)"	"Temperature (K)"	"Speed (ns/day)"	"Time Remaining"
12.5%	1249	-34462.32310300725	296.18992722175085	0	--
25.0%	2498	-34256.932992474976	309.5384422165362	44.6	0:29
37.5%	3747	-34262.00635869912	319.4292600242049	44.2	0:24
50.0%	4996	-34207.6035860422	313.42742950825664	45.5	0:19
INFO [SystemLoader:prepare_simulation] Running Production for 5000 steps, or 10 ps
62.5%	6245	-34094.84048225327	318.3137614730914	42.4	0:15
74.9%	7494	-33912.29320916234	308.6902930492715	39.8	0:10
87.4%	8743	-33681.71865305711	305.08087349583053	37.5	0:05
99.9%	9992	-34039.353253210546	311.1187788403318	36.7	0:00
INFO [SystemLoader:prepare_simulation] Simulation Done. Running MBAR on 50 snapshots.
INFO [SystemLoader:prepare_simulation] Done. Subsampled 30 from 50 snapshots.
INFO [SystemLoader:run_amber] Calculating mmgb/pbsa value...may take awhile.
-25.6433 ± 2.0015 (kcal/mol)
```


#### Citations (to be finished)

[1] Shirts MR and Chodera JD. Statistically optimal analysis of samples from multiple equilibrium states. J. Chem. Phys. 129:124105, 2008 http://dx.doi.org/10.1063/1.2978177

[2] J. D. Chodera, W. C. Swope, J. W. Pitera, C. Seok, and K. A. Dill. Use of the weighted histogram analysis method for the analysis of simulated and parallel tempering simulations. JCTC 3(1):26-41, 2007.

