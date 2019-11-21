# compile all Rmarkdowns
# move the output md file to docs
# remove the html to avoid bloating the repo

# Week 0
cd week0
Rscript -e "rmarkdown::render('R_tutorial.Rmd', output_format = 'github_document')"
mv *.md ../docs/lectures
rm *.html
cd ..

# Week 1
cd week1
Rscript -e "rmarkdown::render('probability_review.Rmd', output_format = 'github_document')"
Rscript -e "rmarkdown::render('basic_visualization.Rmd', output_format = 'github_document')"
mv *.md ../docs/lectures
cp -r *_files ../docs/lectures
rm -r *_files
rm *.html
cd ..

# Week 2
cd week2
Rscript -e "rmarkdown::render('basic_data_wrangling.Rmd', output_format = 'github_document')"
Rscript -e "rmarkdown::render('distributions.Rmd', output_format = 'github_document', envir = new.env())"
mv *.md ../docs/lectures
cp -r *_files ../docs/lectures
rm -r *_files
rm *.html
cd ..

# Week 3
cd week3
Rscript -e "rmarkdown::render('hypothesis_testing.Rmd', output_format = 'github_document')"
Rscript -e "rmarkdown::render('linalg_basics.Rmd', output_format = 'github_document', envir = new.env())"
mv *.md ../docs/lectures
cp -r *_files ../docs/lectures
rm -r *_files
rm *.html
cd ..

# Week 4
cd week4
Rscript -e "rmarkdown::render('linear_models.Rmd', output_format = 'github_document')"
Rscript -e "rmarkdown::render('likelihood.Rmd', output_format = 'github_document', envir = new.env())"
mv *.md ../docs/lectures
cp -r *_files ../docs/lectures
rm -r *_files
rm *.html
cd ..

# Week 5
cd week5
Rscript -e "rmarkdown::render('generalized_linear_models.Rmd', output_format = 'github_document')"
Rscript -e "rmarkdown::render('ANOVA_etc.Rmd', output_format = 'github_document', envir = new.env())"
mv *.md ../docs/lectures
cp -r *_files ../docs/lectures
rm -r *_files
rm *.html
cd ..

# Week 6
cd week6
Rscript -e "rmarkdown::render('model_selection.Rmd', output_format = 'github_document')"
Rscript -e "rmarkdown::render('time_series.Rmd', output_format = 'github_document', envir = new.env())"
mv *.md ../docs/lectures
cp -r *_files ../docs/lectures
rm -r *_files
rm *.html
cd ..

# Week 7
cd week7
Rscript -e "rmarkdown::render('multidimensional_scaling.Rmd', output_format = 'github_document')"
Rscript -e "rmarkdown::render('pca.Rmd', output_format = 'github_document', envir = new.env())"
mv *.md ../docs/lectures
cp -r *_files ../docs/lectures
rm -r *_files
rm *.html
cd ..

# Week 8
cd week8
Rscript -e "rmarkdown::render('clustering.Rmd', output_format = 'github_document')"
Rscript -e "rmarkdown::render('phylo.Rmd', output_format = 'github_document', envir = new.env())"
mv *.md ../docs/lectures
cp -r *_files ../docs/lectures
rm -r *_files
rm *.html
cd ..

