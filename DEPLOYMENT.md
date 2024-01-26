# Deployment Manual

> 📘 **Author:** [Juti Noppornpitak](juti@dnastack.com), [DNAstack](https://dnastack.com)

## Minimum System Requirements

### Server

* Python 3.9
* Python source code for some dependencies
* [Dependencies listed in `requirements.txt`](./requirements.txt)
* 10 GB of RAM

### Deployer's Machine

* Docker
* Terraform
* Access to GCP Project "dnastack-sickkids-strug-lab"
  * Only DNAstack has the access to the GCP project.

## Deployment with Docker and Docker Compose

### Requirements

* A service account with permissions (Storage Object Viewer) to access to the cloud storage bucket

### Getting Started

1. Copy `.env.dist` as `.env`.
   * Here is how to do it:
     ```shell
     cp .env.dist .env
     ```
2. Ensure that the service account's JSON key file exists at `./sa.json`.
2. Run:
   ```shell
   make dev-run
   ```

## Deployment with Cloud Run

1. Run `make image`
2.
