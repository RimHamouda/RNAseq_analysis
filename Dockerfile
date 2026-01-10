FROM rocker/rstudio:4.3.0
# Installation des dépendances système Linux nécessaires pour la bioinfo
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libz-dev \
    libglpk-dev \
    libxt-dev \
    && apt-get clean



