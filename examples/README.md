## Quickstart examples on Summit and Crusher 

### Summit 

- Build SmartRedis following [recipe](https://www.craylabs.org/docs/installation.html#smartsim-on-summit-at-olcf) or download the pre-built [tarbar](https://www.dropbox.com/s/3kxmwj2xepr9tnl/envs.tar?dl=0). TensorFlow is available via open-ce module. 

- Clone the repo
```bash
git clone https://code.ornl.gov/jqyin/deepthermo-wl
```

- Build the binary for DeepThermo 
```bash
cd src
module load gcc
make clean
make -f makefile.summit REDIS_DIR=<path-to-smartredis-built>
cp hea-wl ../examples/summit/
```

- Run the example 
```bash 
cd ../examples/summit 
# edit the env.sh and change REDIS_DIR=<path-to-smartredis-built>
bsub hea-wl.lsf 
```

### Crusher 

- Build TensorFlow using [build-tf-crusher.sh](https://code.ornl.gov/jqyin/deepthermo-wl/-/blob/vae-proposal/utils/build-tf-crusher.sh) or download the pre-built [tarbar](https://www.dropbox.com/s/3kxmwj2xepr9tnl/envs.tar?dl=0).  

- Clone the repo
```bash
git clone https://code.ornl.gov/jqyin/deepthermo-wl
```

- Build the binary for DeepThermo 
```bash
cd src
module load PrgEnv-gnu
make -f makefile.crusher TF_DIR=<path-to-tf-built>
cp hea-wl ../examples/crusher/ 
```

- Run the example 
```bash 
cd ../examples/crusher
# edit the env.sh and change TF_DIR=<path-to-TF-built>
sbatch hea-wl.sb 
```
