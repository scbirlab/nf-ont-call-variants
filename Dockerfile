FROM mambaorg/micromamba:1.5.6

ENV HOME=/home/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/conda
WORKDIR $HOME

USER root
RUN echo "user:x:1001:1001::/home/user:/bin/bash" >> /etc/passwd && \
    mkdir -p /home/user && chown -R 1001:1001 /home/user
RUN mkdir -p $HOME/.conda && chown -R 1000:1000 $HOME
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    gcc \
    g++ \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*
USER 1000

COPY environment.yml /tmp/environment.yml
RUN micromamba create -n env -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH=$MAMBA_ROOT_PREFIX/envs/env/bin:$PATH
