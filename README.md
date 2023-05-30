# MSc-Slurm-Code-Pt-2
This code covers Mafft to ...
If you haven't already, software to install:
Mafft
CIAlign
Optrimal
Julia
IQTREE
raxml-ng
weighted ASTRAL (or ASTRAL III if you do not figure out ASTRAL-w)


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
You will need to use 3 scripts for this and a text file with the cut offs.
OptrimAl, Pasta taster and an R script.

### PASTA Taster
```
#!/bin/bash

while read cutoff_trim
do
        cd /home/DIR/HP_out/alignments/edited
        mkdir $cutoff_trim

        for alignment in *.fasta
        do
          /home/DIR/apps/trimAl/source/trimal -in ${alignment} -out ${cutoff_trim}/${alignment/.fasta}_o.fasta -htmlout ${cutoff_trim}/${alignment/.fasta}_o.html -gt $cutoff_trim

                # check if alignment was trimmed to extinction by trimAl

                if grep ' 0 bp' ${cutoff_trim}/${alignment/.fasta}_o.fasta
                then
                        rm -f ${cutoff_trim}/${alignment/.fasta}_o.fasta
                fi
        done

        cd /home/DIR/HP_out/alignments/edited/${cutoff_trim}
        python3 /home/DIR/apps/AMAS-master/amas/AMAS.py summary -f fasta -d dna -i *_o.fasta

        mv summary.txt ../summary_${cutoff_trim}.txt

done < /home/DIR/edited/cutoff_trim.txt

Rscript /home/DIR/edited/optrimal.R
```

Uploding code
```
scp PASTA_taster.sh DIR/DIR/edited/
```

### OprimAl code
```
#!/bin/bash
#
#SBATCH --chdir=/home/DIR/edited
#SBATCH --job-name=Optrimal
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-user=email@kew.org
#SBATCH --mail-type=END,FAIL

bash PASTA_taster.sh
```

Uploding code
```
scp OptrimAl_script.sh DIR/DIR/edited/
```

Uploding cutoff text
Ex text.
0
0.1
0.3
0.5
0.7
0.9
1
```
scp cutoff_trim.txt DIR/DIR/edited
```
### R script 
Only change pwd
```
# optrimAl.R - run this second

setwd("/home/DIR/HP_out/alignments/edited")
cutoff_trim <- readLines('cutoff_trim.txt')

# create one multiple tables for each threshold value to store AMAS results

amas_table <- read.table('summary_0.txt', header = TRUE)
sites <- data.frame(row.names = amas_table$Alignment_name)
pct <- data.frame(row.names = amas_table$Alignment_name)
filled <- data.frame(row.names = amas_table$Alignment_name)
lost <- data.frame(row.names = amas_table$Alignment_name)

for(i in 1:length(cutoff_trim)){
  amas_table <- read.table(paste('summary_', cutoff_trim[i], '.txt', sep = ''), header = TRUE)
  for(j in amas_table$Alignment_name){
    sites[rownames(sites) == j,i] <- amas_table$Parsimony_informative_sites[amas_table$Alignment_name == j]
    pct[rownames(pct) == j,i] <- as.numeric(amas_table$Proportion_parsimony_informative[amas_table$Alignment_name == j])
    filled[rownames(filled) == j,i] <- amas_table$Total_matrix_cells[amas_table$Alignment_name == j] * (1 - amas_table$Missing_percent[amas_table$Alignment_name == j] / 100)
  }
}

# calculate data loss for each trimming threshold

sites[is.na(sites)] <- 0
pct[is.na(pct)] <- 0

for(i in 1:ncol(filled)){
  lost[,i] <- 1 - filled[,i] / filled[,1]
}

lost[is.na(lost)] <- 1

colnames(sites) <- cutoff_trim
colnames(pct) <- cutoff_trim
colnames(filled) <- cutoff_trim
colnames(lost) <- cutoff_trim

# select optimal trimming threshold
# current criterion is maximum proportion of parsimony informative sites where data loss is no more than one median absolute deviation above the median

optrim <- numeric()
optrim_loss <- numeric()

for(i in rownames(pct)){
  lost_i <- unlist(lost[rownames(lost) == i, ])
  pct_i <- unlist(pct[rownames(pct) == i, ])
  dldp <- data.frame(pct_i, lost_i, row.names = cutoff_trim)
  write.csv(dldp, paste('dldp_', i, '.csv', sep = ''))
  real_loss <- dldp$lost_i[dldp$lost_i < 1]
  diff_loss <- real_loss[2:length(real_loss)] - real_loss[1:(length(real_loss) - 1)]
  median_loss <- median(diff_loss[diff_loss != 0])
  dldp <- subset(dldp, dldp$lost_i <= (median(real_loss) + median_loss))
  if(length(dldp$pct_i) > 0){
    optrim[i] <- rownames(dldp)[dldp$pct_i == max(dldp$pct_i)][[1]]
    optrim_loss[i] <- dldp$lost_i[rownames(dldp) == optrim[i][[1]]]
  } else {
    optrim[i] <- 0
    optrim_loss[i] <- 0
  }
}

# generate graphs to show effect of trimming on informativeness and data loss

for(i in rownames(pct)){
  dldp <- read.csv(paste('dldp_', i, '.csv', sep = ''))
  png(paste('dldp_', i, '.png', sep = ''))
  par(mar = c(5,5,2,5))
  plot(main = i, dldp$lost_i ~ cutoff_trim, ylim = c(0,1), ylab = 'proportion of data lost', xlab = 'strictness of trimming (trimAl gap threshold)', pch = 18, col = 'red')
  par(new = T)
  plot(dldp$pct_i ~ cutoff_trim, xlab = NA, ylab = NA, ylim = c(0,1), axes = F, pch = 16, col = 'blue')
  axis(side = 4)
  mtext(side = 4, line = 3, 'proportion parsimony informative')
  legend(x = 0, y = 1, legend = c('proportion of data lost', 'proportion of parsimony informative sites', 'selected trimming threshold'), pch = c(18, 16, NA), lty = c(NA, NA, 2), col = c('red', 'blue', 'black'), cex = 0.9, bty = 'n')
  if(is.na(optrim[i]) == FALSE){
    lines(c(-0.5, optrim[i]), c(optrim_loss[i], optrim_loss[i]), lty = 2)
    lines(c(-0.5, optrim[i]), c(dldp$pct_i[dldp$X == optrim[i]], dldp$pct_i[dldp$X == optrim[i]]), lty = 2)
    lines(c(optrim[i], optrim[i]), c(-0.5, max(optrim_loss[i], dldp$pct_i[dldp$X == optrim[i]])), lty = 2)
  }
  dev.off()
}

overlost <- names(optrim_loss[optrim_loss > 0.3])

write(overlost, 'overlost.txt', sep = '\n')

file.copy(paste(optrim, '/', names(optrim), sep = ''), getwd())

file.remove(paste(overlost, sep = ''))
```

Uploding R script
```
scp optrimal.R DIR/DIR/HP_out/alignments/edited
```

Running all the scripts by using OptrimAl
```
sbatch /home/DIR/OptrimAl_script.sh
```

If code fails at the R stage run the R script using
```
sbatch /home/gflanaga/scratch/private/test/TMP_Irsalina_Syzygium/HP_out/alignments/edited/optrimal.R

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
## Taper

```
#!/bin/bash
#
#SBATCH --chdir=/home/DIR/HP_out/alignments/edited/alM_r_o
#SBATCH --job-name=taper
#SBATCH --partition=short
#SBATCH --array=1-353%25
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --mail-user=email@kew.org
#SBATCH --mail-type=END,FAIL

echo $SLURM_ARRAY_TASK_ID

name=$(awk -v lineid=$SLURM_ARRAY_TASK_ID 'NR==lineid{print;exit}' /mnt/DIR/HP_out/alignments/test-gene-names-new.txt)

echo $name

/home/DIR/apps/julia-1.9.0/bin/julia /home/DIR/apps/TAPER-master/correction_multi.jl -m N -a N "$name"_alM_r_o.fasta > "$name"_alM_r_cT.fasta
```

















