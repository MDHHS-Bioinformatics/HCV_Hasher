# Use an official Python image with version 3.10.16 or higher
FROM python:3.10.16-slim

# Set the working directory inside the container
WORKDIR /app

# Install ps command and bash
RUN apt-get update && \
    apt-get install -y --no-install-recommends procps bash && \
    rm -rf /var/lib/apt/lists/*

# Copy your script and requirements.txt into the container
COPY hcv_hasher.py /app/
COPY requirements.txt /app/

# Install dependencies from requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
