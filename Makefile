PROD_IMAGE=northamerica-northeast2-docker.pkg.dev/dnastack-sickkids-strug-lab/gwas-svatalog
PROD_VER=$$(git describe | awk -F'[-]' '{print $$1"."$$2}')

.PHONY: image
image:
	@docker buildx create --name ma-builder --bootstrap --use 2> /dev/null || echo "NOTICE: The multiarch builder may have already been configured."
	docker buildx build --push \
		--platform linux/arm64,linux/amd64 \
		--tag $(PROD_IMAGE):$(PROD_VER) .

# For local testing... please copy .env.dist as .env before running this.
.PHONY: dev-run
dev-run:
	docker compose up --build
