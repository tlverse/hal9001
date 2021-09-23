md:
	Rscript -e "rmarkdown::render('README.Rmd', output_file = 'README.md')"

site:
	Rscript -e "rmarkdown::render('README.Rmd', output_file = 'README.md')"
	Rscript -e "pkgdown::build_site()"

check:
	Rscript -e "devtools::check()"

checkfast:
	Rscript -e "devtools::check(build_args = '--no-build-vignettes', run_dont_test = TRUE)"

test:
	Rscript -e "devtools::test()"

doc:
	Rscript -e "devtools::document()"

build:
	Rscript -e "devtools::build()"

buildfast:
	Rscript -e "devtools::build(vignettes = FALSE)"

style:
	Rscript -e "styler::style_pkg()"

pr: style checkfast site
	echo "If all checks have passed, you are ready to submit a new PR."
