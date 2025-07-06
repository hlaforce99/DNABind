# Dockerfile

FROM continuumio/miniconda3

# Install build tools for pip packages (deeptime needs CMake, gcc, make)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
    && rm -rf /var/lib/apt/lists/*

COPY environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml

SHELL ["conda", "run", "-n", "dnarligand-freeenergy", "/bin/bash", "-c"]

WORKDIR /app
COPY . /app

CMD ["bash"]

