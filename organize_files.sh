# compile all Rmarkdowns
# move the output md file to docs
# remove the html to avoid bloating the repo

# Week 0
cd week0
Rscript -e "rmarkdown::render('R_tutorial.Rmd')"
mv *.md ../docs/lectures
rm *.html
cd ..

# Week 1
cd week1
Rscript -e "rmarkdown::render('probability_review.Rmd')"
Rscript -e "rmarkdown::render('basic_visualization.Rmd')"
mv *.md ../docs/lectures
cp -r *_files ../docs/lectures
rm -r *_files
rm *.html
cd ..

# Week 2
cd week1
Rscript -e "rmarkdown::render('basic_data_wrangling.Rmd')"
#Rscript -e "rmarkdown::render('basic_visualization.Rmd')"
mv *.md ../docs/lectures
cp -r *_files ../docs/lectures
rm -r *_files
rm *.html
cd ..


