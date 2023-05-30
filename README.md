# MSc-Slurm-Code-Pt-2
This code covers Mafft to ...

## MAFFT
```
#!/bin/bash
#
#SBATCH --chdir=/home/DIR/HP_out
#SBATCH --job-name=mafft
#SBATCH --partition=medium
#SBATCH --array=1-353      
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --mail-user=*email@kew.org*
#SBATCH --mail-type=END,FAIL


echo $SLURM_ARRAY_TASK_ID

name=$(awk -v lineid=$SLURM_ARRAY_TASK_ID 'NR==lineid{print;exit}' /mnt/DIR/HP_out/test-gene-names.txt)

echo $name

mafft --thread 8 --genafpair --adjustdirectionaccurately --maxiterate 1000 /home/DIR/HP_out/genes/"$name".FNA > "$name"_alM.fasta
```
-----------------
Uploding code
```
scp mafft-script.sh DIR/DIR/apps
```

Running script
```
sbatch /home/DIR/apps/mafft-script.sh
```
Remove any "_R_" once the program is complete
```
for f in *_alM.fasta; do (sed -e 's/_R_//g' $f > ${f/.fasta}_r.fasta); done
```
## Optrimal

Uploding code
```
scp PASTA_taster.sh DIR/DIR/edited/
```

Uploding cutoff text
```
scp cutoff_trim.txt DIR/DIR/edited
```

Uploding code
```
scp OptrimAl_script.sh DIR/DIR/edited/
```

Running script
```
sbatch /home/DIR/OptrimAl_script.sh
```

----------------------
## CIAlign
```
#!/bin/bash
#
#SBATCH --chdir=/home/DIR/alM_r_o
#SBATCH --job-name=CIAlign
#SBATCH --partition=short
#SBATCH --array=1-353 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --mail-user=*email*@kew.org
#SBATCH --mail-type=END,FAIL

export PATH=$PATH:/home/gflanaga/scratch/apps/CIAlign-master/CIAlign/

echo $SLURM_ARRAY_TASK_ID

name=$(awk -v lineid=$SLURM_ARRAY_TASK_ID 'NR==lineid{print;exit}' /mnt/DIR/HP_out/test-gene-names.txt)

echo $name

CIAlign.py --infile "$name"_alM_r_o.fasta --outfile_stem "$name"_alM_r_o_CI85 --remove_divergent --remove_divergent_minperc 0.85 --retain_str Ceroxylon --plot_input --plot_output --plot_markup
mv "$name"_alM_r_o_CI85_cleaned.fasta "$name"_alM_r_o_CI85.fasta
```

Uploading Code
```
scp CIAlign_slurm.sh DIR/apps
```

Run
```
sbatch /home/DIR/apps/CIAlign_slurm.sh 
```
----------------------------------------------------------------------


















