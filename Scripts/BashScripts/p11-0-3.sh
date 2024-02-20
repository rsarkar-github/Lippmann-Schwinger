# Run from within this directory

cd /home/rsarkar/
source .bashrc
conda activate py39
cd /home/rsarkar/Research/Thesis
python -m Lippmann-Schwinger.Scripts.p11_initial_vel_solves 0 3 0
python -m Lippmann-Schwinger.Scripts.p11_initial_vel_solves 0 3 1
python -m Lippmann-Schwinger.Scripts.p11_initial_vel_solves 0 3 2
