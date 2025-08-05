FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

# Install generic dependencies
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y --fix-missing \
    software-properties-common \
    build-essential \
    libxrender1 \
    libxtst6 \
    ncbi-blast+ \
    pymol \
    wget && \
    rm -rf /var/lib/apt/lists/*

# Install Github CLI
RUN (type -p wget >/dev/null || (apt update && apt install wget -y)) \
	&& mkdir -p -m 755 /etc/apt/keyrings \
	&& out=$(mktemp) && wget -nv -O$out https://cli.github.com/packages/githubcli-archive-keyring.gpg \
	&& cat $out | tee /etc/apt/keyrings/githubcli-archive-keyring.gpg > /dev/null \
	&& chmod go+r /etc/apt/keyrings/githubcli-archive-keyring.gpg \
	&& mkdir -p -m 755 /etc/apt/sources.list.d \
	&& echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
	&& apt update \
	&& apt install gh -y

# Install Python and Hatch for package management
RUN add-apt-repository ppa:deadsnakes/ppa && apt-get update && apt-get install -y python3.12 python3.12-dev python3-pip
RUN pip install --upgrade pip
RUN pip install hatch
ENV PYTHONUNBUFFERED=1

# Set up data and variables
ENV INSTALL_PATH=/mite_web
RUN mkdir -p $INSTALL_PATH
WORKDIR $INSTALL_PATH
COPY ./mite_web/pyproject.toml .
RUN mkdir -p mite_web && touch mite_web/__init__.py

# Create the environment, download Python dependencies
RUN hatch env create

# Copy the rest of the files
COPY ./mite_web .

# Populate the app with data
RUN hatch run python3.12 mite_web/prepare_mite_data.py

# Create the enzyme images
RUN chmod +x ./mite_web/run_pymol.sh
RUN ./mite_web/run_pymol.sh

EXPOSE 8004

RUN chmod +x ./entrypoint_docker.sh
ENTRYPOINT ["./entrypoint_docker.sh"]
