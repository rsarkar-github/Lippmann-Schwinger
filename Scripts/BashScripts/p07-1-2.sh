# Run from within this directory

cd /home/rsarkar/
source .bashrc
conda activate py39
cd /home/rsarkar/Research/Thesis
python -m Lippmann-Schwinger.Scripts.p07_lse_solves 1 2 0
python -m Lippmann-Schwinger.Scripts.p07_lse_solves 1 2 1
python -m Lippmann-Schwinger.Scripts.p07_lse_solves 1 2 2
