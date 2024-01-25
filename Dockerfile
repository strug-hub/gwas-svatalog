##### Approach: Python #####
# FROM python:3.11

##### Approach: Debian #####
FROM debian:bookworm

ENV GCSFUSE_REPO gcsfuse-bookworm
ENV GOOGLE_REPO_KEY_FILE "/usr/share/keyrings/cloud.google.asc"

RUN apt update \
    && apt install -qy \
        apt-transport-https \
        ca-certificates \
        curl \
        gnupg \
        python3 \
        python3-dev \
        python3-pip
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | tee ${GOOGLE_REPO_KEY_FILE}
RUN echo "deb [signed-by=${GOOGLE_REPO_KEY_FILE}] https://packages.cloud.google.com/apt cloud-sdk main" \
    | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
# RUN echo "deb [signed-by=${GOOGLE_REPO_KEY_FILE}] https://packages.cloud.google.com/apt $GCSFUSE_REPO main" \
#     | tee /etc/apt/sources.list.d/gcsfuse.list
RUN apt update \
    && apt install -qy \
        google-cloud-cli

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
