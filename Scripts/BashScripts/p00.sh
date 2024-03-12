# Run from within this directory

cd ~
source .bashrc
conda activate py39

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo $SCRIPT_DIR
cd $SCRIPT_DIR ..
#cd /home/rsarkar/Research/Thesis
python -m Lippmann-Schwinger.Scripts.p01_create_vz_marmousi