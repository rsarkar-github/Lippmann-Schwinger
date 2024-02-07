# Run from within this directory

cd /home/rsarkar/
source .bashrc
conda activate py39
cd /home/rsarkar/Research/Thesis
python -m Lippmann-Schwinger.Scripts.p05_calculate_green_func 0 0
python -m Lippmann-Schwinger.Scripts.p05_calculate_green_func 0 1
python -m Lippmann-Schwinger.Scripts.p05_calculate_green_func 0 2
python -m Lippmann-Schwinger.Scripts.p05_calculate_green_func 1 0
python -m Lippmann-Schwinger.Scripts.p05_calculate_green_func 1 1
python -m Lippmann-Schwinger.Scripts.p05_calculate_green_func 1 2
python -m Lippmann-Schwinger.Scripts.p05_calculate_green_func 2 0
python -m Lippmann-Schwinger.Scripts.p05_calculate_green_func 2 1
python -m Lippmann-Schwinger.Scripts.p05_calculate_green_func 2 2