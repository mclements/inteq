R = R-devel

test:
	$(R) --slave -e "devtools::test()"

check: build
	$(R) CMD check --as-cran inteq_`awk '/Version:/ {print $$2}' DESCRIPTION`.tar.gz

build:
	$(R) CMD build .
