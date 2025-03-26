# Use an official Python image with version 3.9.19 or higher
FROM python:3.9.19-slim

# Set the working directory inside the container
WORKDIR /app

# Copy the Python script and any other necessary files into the container
COPY hcv_hasher.py /app/

# Install required dependencies
RUN pip install --no-cache-dir \
    mmh3 \
    pandas \
    biopython \
    numpy

# Set the default command to run the script
ENTRYPOINT ["python", "hcv_hasher.py"]
