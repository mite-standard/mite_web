FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

# Get required libraries
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y \
    software-properties-common \
    build-essential \
    && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get upgrade -y
RUN apt-get install -y libxrender1
RUN apt-get install -y libxtst6

# Get NCBI-BLAST+
RUN apt-get install -y ncbi-blast+

# Get PyMol
RUN apt-get install -y pymol

# Install Python and Hatch for package management
RUN add-apt-repository ppa:deadsnakes/ppa && apt-get update && apt-get install -y python3.12 python3.12-dev python3-pip
RUN pip install --upgrade pip
RUN pip install hatch
ENV PYTHONUNBUFFERED 1

# Set up data and variables
ENV INSTALL_PATH /mite_web
RUN mkdir -p $INSTALL_PATH
WORKDIR $INSTALL_PATH
COPY ./mite_web .

# Create the environment, download Python dependencies
RUN hatch env create

# Populate the app with data
RUN hatch run python3.12 mite_web/prepare_mite_data.py

# Create the enzyme images
RUN chmod +x ./mite_web/run_pymol.sh
RUN ./mite_web/run_pymol.sh

RUN chmod +x ./entrypoint_docker.sh

EXPOSE 8004

CMD ["./entrypoint_docker.sh"]
