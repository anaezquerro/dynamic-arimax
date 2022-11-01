FROM ubuntu:latest
ARG DEBIAN_FRONTEND=noninteractive


WORKDIR /app
COPY . /app/

RUN apt-get -y update && apt-get -y upgrade

RUN apt -y install wget gnupg make python3 python3-pip && apt -y update && apt -y upgrade && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/r-project.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" | tee -a /etc/apt/sources.list.d/r-project.list && \
    apt -y update && apt -y install --no-install-recommends r-base r-base-dev r-cran-devtools

RUN Rscript installation.R
RUN Rscript test.R

