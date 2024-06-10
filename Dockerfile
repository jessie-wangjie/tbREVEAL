# Start with the micromamba base image
FROM mambaorg/micromamba:1.5.1

USER root

# Set the maintainer label
LABEL maintainer="thomas.biondi@tome.bio"

# Set the working directory
WORKDIR /app

# Copy your environment file to the docker image
COPY environment.yml .

# Install necessary packages using apt-get
RUN apt-get update && apt-get install -y wget bzip2 procps python3.10 unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Use micromamba to create an environment from the environment.yml file
RUN micromamba install -y -n base -f environment.yml && micromamba clean --all --yes

COPY bin/ ./bin/
RUN chmod 755 ./bin/*
COPY root.crt  /root/.postgresql/root.crt

# Make RUN commands use the base environment
SHELL ["micromamba", "run", "-n", "base", "/bin/bash", "-c"]

# Install AWS CLI
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm -rf awscliv2.zip aws

# Install BaseSpace CLI
RUN wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs" -O /usr/local/bin/bs && \
    chmod u+x /usr/local/bin/bs

# Use this if you want to enter the container at the base environment
ENTRYPOINT ["micromamba", "run", "-n", "base"]
