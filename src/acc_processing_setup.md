ssh mulqueen@acc.ohsu.edu
ssh arc-infra-1
screen
srun --partition=guest --cpus-per-task=30 --time=12:00:00 --mem=300G --nodes=1 --pty /bin/bash
#or srun --partition=interactive --cpus-per-task=30 --time=12:00:00 --mem=300G --nodes=1 --pty /bin/bash
#or srun --partition=cedar --cpus-per-task=30 --time=12:00:00 --mem=300G --nodes=1 --pty /bin/bash

sif="/home/groups/CEDAR/mulqueen/bc_multiome/multiome_bc.sif"
singularity shell --bind /home/groups/CEDAR/mulqueen/bc_multiome --bind /home/groups/MohammedLab --bind /home/groups/CEDAR $sif

#tunnel for sftp
ssh -L 2222:arc-infra-1:22 mulqueen@acc.ohsu.edu
#in new terminal
sftp -P 2222 mulqueen@localhost