##### Approach: Python #####
# FROM python:3.11

##### Approach: Debian #####
FROM debian:bookworm

RUN apt update \
    && apt install -qy \
        apt-transport-https \
        ca-certificates \
        curl \
        gnupg \
        python3 \
        python3-dev \
        python3-pip

##### Approach: ALL #####

EXPOSE 1234

RUN useradd --home-dir /app -m --user-group -u 2000 app

WORKDIR /app

RUN chown app /app

# Install dependencies
COPY requirements.txt .
RUN pip3 install --break-system-packages --no-warn-script-location -q -r requirements.txt
RUN pip3 install --break-system-packages --no-warn-script-location -q crcmod  # To speed up gsutil operations.
RUN rm requirements.txt

# Copy the other files
ADD startup .
ADD LICENSE .
ADD assets ./assets
ADD app.py .

##### Approach: Python #####
# CMD python3 app.py

USER app

##### Approach: Debian #####
CMD ./startup
