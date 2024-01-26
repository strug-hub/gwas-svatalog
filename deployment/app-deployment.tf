variable "deployment_project" {
  description = "Deployment Project"
  default     = "dnastack-sickkids-strug-lab"
}

variable "deployment_region" {
  description = "Deployment Region"
  default     = "northamerica-northeast2"
}

variable "app_version" {
  description = "App Version"
  default     = "1.0.1"
}

variable "lb_ssl" {
  description = "Load Balancer SSL"
  default     = true
}

variable "lb_name" {
  description = "Load Balancer Name"
  default     = "svatalog-lb"
}

variable "cloud_dns_zone" {
  description = "DNS zone to link to Cloud Run service"
  default     = "strug-lab-svatalog"
}

variable "cloud_dns_entry" {
  description = "DNS entry to link to Cloud Run service"
  default     = "svatalog.research.sickkids.ca"
}

terraform {
  backend "gcs" {
    bucket = "sickkids-gwas-svatalog-terraform"
    prefix = "terraform/state"
  }
  required_providers {
    # DOCS: https://registry.terraform.io/providers/hashicorp/google/latest
    google = {
      source  = "hashicorp/google"
      version = "4.51.0"
    }
  }
}

provider "google" {
  project = var.deployment_project
  region  = var.deployment_region
}

module "lb-http" {
  # DOCS: https://github.com/terraform-google-modules/terraform-google-lb-http/tree/master/modules/serverless_negs
  source  = "terraform-google-modules/lb-http/google//modules/serverless_negs"
  version = "10.0.0"

  name    = var.lb_name
  project = var.deployment_project

  ssl                             = var.lb_ssl
  managed_ssl_certificate_domains = [var.cloud_dns_entry]
  https_redirect                  = var.lb_ssl
  labels                          = { "maintainer" = "dnastack", "platform" = "cloud-run" }

  backends = {
    default = {
      description = null

      groups = [
        {
          group = google_compute_region_network_endpoint_group.app_neg.id
        }
      ]
      enable_cdn = false

      iap_config = {
        enable = false
      }

      log_config = {
        enable = false
      }
    }
  }
}

resource "google_compute_region_network_endpoint_group" "app_neg" {
  project               = var.deployment_project
  provider              = google-beta
  name                  = "svatalog-neg"
  network_endpoint_type = "SERVERLESS"
  region                = var.deployment_region
  cloud_run {
    service = google_cloud_run_v2_service.app.name
  }
}

resource "google_cloud_run_v2_service" "app" {
  project  = var.deployment_project
  location = var.deployment_region
  name     = "sickkids-gwas-svatalog"
  client   = "terraform"

  template {
    scaling {
      min_instance_count = 1
      max_instance_count = 10
    }
    containers {
      image = "northamerica-northeast2-docker.pkg.dev/dnastack-sickkids-strug-lab/webapp/gwas-svatalog:${var.app_version}"

      resources {
        limits = {
          cpu    = 4
          memory = "11Gi"
        }
      }

      # Enable the DNAstack debug mode to monitor memory usage (for now)
      env {
        name  = "DNASTACK_DEBUG"
        value = "true"
      }

      # Disable the Dash/Flask debug mode to also disable the auto-reload feature as the data won't be updated live.
      env {
        name  = "SVATALOG_DEBUG_DISABLED"
        value = "true"
      }

      startup_probe {
        initial_delay_seconds = 180 # 3 minutes delay
        timeout_seconds       = 60
        period_seconds        = 15
        failure_threshold     = 5
        http_get {
          path = "/"
        }
      }

      liveness_probe {
        initial_delay_seconds = 300 # 5 minutes delay
        timeout_seconds       = 30
        period_seconds        = 15
        failure_threshold     = 5
        http_get {
          path = "/"
        }
      }
    }
  }
}

resource "google_cloud_run_v2_service_iam_member" "noauth" {
  project  = google_cloud_run_v2_service.app.project
  location = google_cloud_run_v2_service.app.location
  name     = google_cloud_run_v2_service.app.name
  role     = "roles/run.invoker"
  member   = "allUsers"
}

output "cloud_run_default_app_url" {
  value = google_cloud_run_v2_service.app.uri
}

resource "google_dns_record_set" "app-dns-caa" {
  name         = "${var.cloud_dns_entry}."
  type         = "CAA"
  ttl          = 3600
  managed_zone = var.cloud_dns_zone

  rrdatas = ["0 issue \"letsencrypt.org\""]
}

resource "google_dns_record_set" "app-dns" {
  name         = "${var.cloud_dns_entry}."
  type         = "A"
  ttl          = 3600
  managed_zone = var.cloud_dns_zone

  rrdatas = [module.lb-http.external_ip]
}
