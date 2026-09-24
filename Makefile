IMAGE ?= dna-compare
PORT ?= 8765
DATA_DIR ?= $(CURDIR)/data

.PHONY: build serve

build:
	docker build -t $(IMAGE) .

serve: build
	docker run --rm --name $(IMAGE) -p $(PORT):8765 \
		-v $(DATA_DIR):/app/data \
		$(IMAGE)
