FROM nvidia/cuda:11.1-base-ubuntu20.04
RUN apt update && apt install -y python3.8 && apt install -y python3-pip
RUN apt install -y vim
RUN apt install -y git
RUN mkdir /configfiles
COPY requirements.txt /configfiles
RUN pip3 install -r /configfiles/requirements.txt
RUN pip3 install torch==1.10.1+cu113 -f https://download.pytorch.org/whl/cu113/torch_stable.html
