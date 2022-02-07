
# beataml2_manuscript


git clone https://github.com/biodev/beataml2_manuscript.git beataml2_manuscript

## get data from xxxx

Should contain:

1-s2.0-S0092867419300947-mmc3.xlsx
beataml_waves1to4_norm_exp_dbgap.txt	
beataml_probit_curve_fits_v4_dbgap.txt
beataml_wes_wv1to4_mutations_dbgap.txt
Table_S4.xlsx
Table_S1.xlsx
wgcna/merged_older_wgcna_kme.RData			

cp -r xxx/data  beataml2_manuscript/

## To run

```
cd beataml2_manuscript

R
renv::restore()
#note it mentions repair of dependency tree, seems to be an renv issue with packages not directly used
targets::tar_make()

```

### Note: To reproduce fully, need to use R-4.03

To isolate this procedure from your system's R, use of Docker is recommended.

## Docker instructions

Use of the `no_font` branch is recommended as the specification of Arial font family can cause issues in Ubuntu

```
cd beataml2_manuscript
git checkout no_font
cd ..

docker pull r-base:4.0.3

docker run -it -w /data/beataml2_manuscript -v $PWD:/data r-base:4.0.3 /bin/bash

apt-get -y update
apt-get -y  install libcairo2-dev libxt-dev libxml2-dev libssl-dev libcurl4-openssl-dev

mkdir figures
mkdir output_files

```

$ R
```
renv::restore()
targets::tar_make()
```