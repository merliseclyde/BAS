usethis::use_release_issue()

# check email for results
devtools::check_win_devel()
devtools::check_win_release()

# mac_builder M1 mac 
# https://mac.r-project.org/macbuilder/submit.html

# checks for rhub

ch <- rhub::check_for_cran(".", show_status = FALSE)
ch$cran_summary()
ch$update() 

# check for M1Mac
url('https://mac.r-project.org/macbuilder/submit.html')

# to submit to CRAN

devtools::submit_cran()
