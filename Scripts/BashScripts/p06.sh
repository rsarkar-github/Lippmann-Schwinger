# Run from within this directory

cd /home/rsarkar/
source .bashrc
conda activate py39
cd /home/rsarkar/Research/Thesis
python -m Lippmann-Schwinger.Scripts.p06_create_sources 0
python -m Lippmann-Schwinger.Scripts.p06_create_sources 1
python -m Lippmann-Schwinger.Scripts.p06_create_sources 2
