# Run from within this directory

cd /home/rsarkar/
source .bashrc
conda activate py39
cd /home/rsarkar/Research/Thesis
python -m Lippmann-Schwinger.Scripts.p08_helmholtz_solves 2 2 0
python -m Lippmann-Schwinger.Scripts.p08_helmholtz_solves 2 2 1
python -m Lippmann-Schwinger.Scripts.p08_helmholtz_solves 2 2 2
