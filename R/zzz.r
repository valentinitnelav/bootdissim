# Declare global variables and native symbol objects ----------------------

# Doing so, avoids the note from devtools::check():
# "...no visible global function definition for...". or
# "...no visible binding for global variable..."
# See https://stackoverflow.com/a/12429344/5193830
# or https://stackoverflow.com/a/17807914/5193830

.onLoad <- function(...) {
  if (getRversion() >= "4.0")
    utils::globalVariables(
      c(
        'spl_size',
        'simulation_id',
        'value',
        'ci_low',
        'ci_up',
        'counts',
        '.N'
      )
    )
}


# Package startup message -------------------------------------------------

.onAttach <- function(...) {
  packageStartupMessage(strwrap("Bootstrap dissimilarity metrics between ecological networks",
                                indent = 5))
}
