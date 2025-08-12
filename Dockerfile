# Use an official Python image with version 3.9.19 or higher
FROM python:3.9.19-slim

# Set the working directory inside the container
WORKDIR /app

# Install ps command and bash
RUN apt-get update && \
    apt-get install -y procps bash && \
    rm -rf /var/lib/apt/lists/*

COPY hcv_hasher.py /app/

# Install required dependencies
RUN pip install --no-cache-dir \
    mmh3 \
    pandas \
    biopython \
    numpy


