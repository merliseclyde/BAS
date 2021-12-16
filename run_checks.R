# check email fro results
devtools::check_win_devel()
devtools::check_win_release()

# checks for rhub

ch <- rhub::check_for_cran(".", show_status = FALSE)
ch$cran_summary()
ch$update() 
