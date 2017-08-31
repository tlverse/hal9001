md:
	Rscript -e "rmarkdown::render('README.Rmd')"

site:
	Rscript -e "pkgdown::build_site()"

check:
	Rscript -e "devtools::check()"

test:
	Rscript -e "devtools::test()"

doc:
	Rscript -e "devtools::document()"

