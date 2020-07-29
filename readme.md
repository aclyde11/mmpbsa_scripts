 
 
 
 ### Setup
 
 Please download and install Amber20. At the end, you should have a path to amber and set the environment variable:
 ```
export AMBERHOME=/Users/austin/Downloads/amber20
```
where inside the folder there is a file ```amber.sh```. Do not source it in the terminal where you will run the code, however.


 
Create a new Conda env with python 3.7 and install the following packages below.
```
conda create -n myenv python=3.7 -y
conda activate myenv
conda install -c conda-forge -c omnia -c openeye openmm openeye-toolkits
export OE_LICENSE=...
``` 


### Usage
There will be more features added in the future. For now, this is the code on the master branch. 

Please make sure your working directory is the base of this repo.

```
cd path-to-mmgbsa_scripts/
```

Here is an example. This will only output the final GBSA value as v=0, v of 1 or 2 will output more information. This will relax
the system for 100,000 steps (step size is 2fs). Following that, there is 100ps of a production run, where 5 snapshots are taken.
These five snapshots are used to compute delta g values, and averaged over those five. 

PBSA is very slow, don't go over five...
```
CUDA_VISIBLE_DEVICES=0 python main.py --amber_path /Users/austin/Downloads/amber20 --pdb /Users/austin/MPro-docked/MPro_0387_Gen3L_5.pdb \ 
               --odir testdir -v 0 --ps 100 --traj_frames 5 --method gbsa --equil_steps 100000 --use_gpu
```