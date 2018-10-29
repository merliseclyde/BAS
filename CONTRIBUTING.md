# Contributing to BAS development

The goal of this guide is to help you contribute to BAS. The guide is divided into three main pieces:

1. Filing a [bug report](https://github.com/merliseclyde/BAS/issues/new?template=bug_report.md) or [feature request](https://github.com/merliseclyde/BAS/issues/new?template=feature_request.md) in an issue on Github.

2. Suggesting a change via a pull request.

3. Coding Style Guide for Contributions to BAS

## Issues

### Feature Requests

Do you wish that you could easily extract additional information from BAS objects ?  Or would you like to see new functionality available in BAS?   If so, feel free to fill out a [feature request](https://github.com/merliseclyde/BAS/issues/new?template=feature_request.md)! Please describe in as much detail what you would like added.  This can be  anything from just the idea on up to code if you are a  more advanced user!    

### Bug Reports

When filing a [bug report](https://github.com/merliseclyde/BAS/issues/new?template=bug_report.md) for an issue, the most important thing is to include a minimal reproducible example so that I can quickly verify the problem, and then figure out how to fix it. There are three things you need to include to make your example reproducible: required packages, data, code.

1.  **Packages** should be loaded at the top of the script, 
so it's easy to see which ones the example needs.

2.  The easiest way to include **data** is to use `dput()` to generate the R code
    to recreate it. For example, to recreate the `mtcars` dataset in R,
    I'd perform the following steps:

       a. Run `dput(mtcars)` in R
       b. Copy the output
       c. In my reproducible script, type `mtcars <- ` then paste.

    But even better is if you can create a `data.frame()` with just a handful
    of rows and columns that still illustrates the problem.

3.  Spend a little bit of time ensuring that your **code** is easy for others to
    read:

    * make sure you've used spaces and your variable names are concise,
    but informative (I am OK with using ".", camel case or "_" in variable names to improve
    readibility.  For more details see the [Style Guide](STYLE_GUIDE.html)

    * use comments to indicate where your problem lies

    * do your best to remove everything that is not related to the problem.
     The shorter your code is, the easier it is to understand.

You can check you have actually made a reproducible example by starting up a fresh R session and pasting your script in.

(Unless you've been specifically asked for it, please don't include the output of `sessionInfo()`.)

### Other issues

If you are not sure something is a bug or an undocumented feature, see possible errors in the help files or documentation that could use clarification (or any other issue) please file a regular issue

## Pull requests

To contribute a change to `BAS`, follow these steps:

1. Create a branch in git and make your changes, ideally using `commit -s` to sign-off on commits  under the [Developer Certificate of Origin](https://developercertificate.org).  Make sure that your branch passes `R CMD check`. 

2. Push branch to github and issue pull request (PR).

3. Discuss the pull request.

4. Iterate until either we accept the PR or decide that it's not
   a good fit for `BAS`.

Each of these steps are described in more detail below. This might feel overwhelming the first time you get set up, but it gets easier with practice. If you get stuck at any point, please reach out for help.

If you're not familiar with git or github, please start by reading <http://r-pkgs.had.co.nz/git.html>


Pull requests will be evaluated against the following checklist:

1.  __Motivation__. Your pull request should clearly and concisely
motivates the need for change. Please describe the problem and show
how your pull request solves it as concisely as possible.

Also include this motivation in `NEWS` so that when a new release of
`BAS` comes out it's easy for users to see what's changed. Add your
item at the top of the file and use markdown for formatting. The
news item should end with `(@yourGithubUsername, #the_issue_number)`.

2.  __Only related changes__. Before you submit your pull request,
please check to make sure that you haven't accidentally included any
unrelated  changes. These make it harder to see exactly what's changed,
and to evaluate any unexpected side effects.

  Each PR corresponds to a git branch, so if you expect to submit
  multiple changes make sure to create multiple branches. If you have
  multiple changes that depend on each other, start with the first one
  and don't submit any others until the first one has been processed.


3.  __Document__ If you're adding new parameters or a new function,
you'll also need to document them with [roxygen](https://github.com/klutometis/roxygen).   Please add a short
example to the appropriate function and optionally in the package vignettes. Make sure to re-run `devtools::document()` on the code before submitting.  (and be sure to include your name in the authors for the function!)


4.  __Testing__ If fixing a bug or adding a new feature, you should add a [testthat](https://github.com/hadley/testthat) unit test.



This seems like a lot of work but don't worry if your pull request isn't perfect. It's a learning process and I will be on hand to help you out. A pull request is a process, and unless you've submitted a few in the past it's unlikely that your pull request will be accepted as is.

Please don't submit pull requests that change existing behaviour. Instead, think about how you can add a new feature in a minimally invasive way.


## Style Guild for Contributing to BAS

A consistent style improves the readibility of your code.
I am not wed to a particular style and generally draw from both 
[Google Style Guide](https://google.github.io/styleguide/Rguide.xml) as well as
[Hadley Wickham's Style Guide](http://adv-r.had.co.nz/Style.html) as noted 
below.  Using a package such as [styler](http://styler.r-lib.org)  to 
enforce styling based on the TidyVerse is helpful, but is not required.

1. *Function and Variable names*: Use informative  names using ".", camel case,
or "_" to improve readibility, i.e. `variable.name`, `VariableName` or
`variable_name`, rather than `foo` or `xxx`.  I tend to avoid `_` for historical
reasons going back to `S`.  

2. *Assignment*: Use either  `<-`  or `=` for assignment but be consistent
within your contribution.  As many style guides prefer `<-` I suggest using
`styler` to enforce use of `<-`, but I
am OK with `=` as it is shorter to type!  Just be consistent within the contributed code.

3. *Spaces*: Include spaces around operators `=, +, -, <-,  ==` etc to improve
readibility.  Put a space after a comman, but not before it.   `:` or `::` never
have spaces around them.  Additional spaces or newlines are fine if they
improve readibilty of code (e.g. aligmnent of arguments).

4. *Comments*: Comment your code whenever you can. Explain the why and not the
what as that should be clear from your code. Use `#` to start a comment followed
by a space and capitalize the first letter; short inline comments (comments on 
the same line as code) need two spaces before the `#`

5.  *Curly Braces*:  An opening curly brace should never go on its own line,
while a closing curly brace should go on its own line.  An exception is with a
short conditional statement shuch as `else` where the code may fit on one line.

6.  *Indentation*:  Use 2 spaces rather than tabs per level of indentation. 
Indent code inside  curly braces.

7.  *Semi-colons*: Do not use semi-colons to put more than one statement on a
line.

8. *Line Length*:  Use up to 80 characters per line of code. RStudio has a 
setting to display a vertical line at 80 characters to visually assist you.
Turn on by going to

  Tools -> Global Options… -> Code -> Display -> Show margin.
  
9.  File names for `R` code should be informative and end in `.R`.  Use `-` to improve readibility. Do not include any spaces in file names!   
  


_Contributing was adopted from ggplot2's_
[CONTRIBUTING.md](https://github.com/tidyverse/ggplot2/blob/master/CONTRIBUTING.md)
