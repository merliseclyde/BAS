# first merge main devel in terminal
# on devel
git merge main
# change to main
git checkout main
# merge changes from devel into main
git merge devel


# Check current CRAN check results
usethis::use_release_issue()


urlchecker::url_check()
devtools::build_readme()
devtools::check(remote = TRUE, manual = TRUE)

devtools::install_github("r-lib/revdepcheck")

revdepcheck::revdep_reset()
revdepcheck::revdep_check(num_workers = 4)
revdepcheck::revdep_report_cran()

# check email for results
devtools::check_win_devel()


# check rhub. (see github actions to trigger rhub workflow
# rhub2::rhub_check(platforms = "valgrind")
# rhub2::rhub_check(platforms = "valgrind", branch="devel")



# to submit to CRAN
usethis::use_version('patch')

devtools::submit_cran()

usethis::use_dev_version(push = TRUE)
