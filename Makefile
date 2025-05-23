.PHONY: run test

TEST_EXEC=build/main_test

run:
	chmod +x run.sh && ./run.sh

test: build
	@echo "=== Running tests ==="
	$(TEST_EXEC)

clean:
	rm -rf build