#!/bin/bash
#Bash strartup
#installs nest,nestml,pip3
sudo add-apt-repository ppa:nest-simulator/nest
sudo apt-get update
sudo apt-get install nest
sudo apt-get install libssl-dev
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null |
    sudo apt-key add -
sudo apt-add-repository 'deb https://apt.kitware.com/ubuntu/ bionic main'
sudo apt update
sudo apt install cmake
#sudo apt install cmake ubnutu has an old cmake
sudo apt install python3-pip
alias python=python3
alias pip=pip3
pip install numpy,scipy,ipython, seaborn,networkx,h5py,cython,html
cd /home/$USER
git clone  https://github.com/nest/nestml.git
cd /home/$USER/nestml
sudo python3 -m pip install -r requirements.txt
sudo python3 setup.py install
#echo "export LD_LIBRARY_PATH=/home/$USER/.local/lib/" >> ~/.bashrc
pip install pygsl
