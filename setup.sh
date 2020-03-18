#!/bin/bash

sudo add-apt-repository ppa:nest-simulator/nest
sudo apt-get update
sudo apt-get install nest
cd /home/$USER
git clone  https://github.com/nest/nestml.git
cd /home/$USER/nestml
sudo python3 -m pip install -r requirements.txt
sudo python3 setup.py install
echo "export LD_LIBRARY_PATH=/home/$USER/.local/lib/" >> ~/.bashrc
alias python=python3
