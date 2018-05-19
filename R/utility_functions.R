#' Checks whether a variable exists and is not NULL
#'
#' @param varname variable name
#'
#' @return TRUE or FALSE
#' @export
enn = function(varname) {
  if(exists(as.character(substitute(varname)))) if(!is.null(varname)) return(TRUE)
  return(FALSE)
}


## Function to add rows to a data frame
insertRow2 = function(existingDF, newrow, rownum, newname="new_row", ...) {
  existingDF = rbind(existingDF, newrow)
  rownames(existingDF)[NROW(existingDF)] = newname
  existingDF = existingDF[order(c(1:(nrow(existingDF)-1),rownum-0.5)),]
  return(existingDF)
}


#' Function to trim leading and trailing whitespace
#'
#' @param x string to be trimmed
#'
#' @return Trimmed string
#' @export
trimwhite = function (x) gsub("^\\s+|\\s+$", "", x)


#' Fail-tolerant function to convert a string to numeric if applicable
#'
#' @param string string with potential numeric values
#'
#' @return Numeric variable (if appropriate), otherwise string
#' @export
string2num = function(string) {
  if (suppressWarnings(all(!is.na(as.numeric(as.character(string)))))) {
    return(as.numeric(as.character(string)))
  } else {
    return(string)
  }
}


## Checks if all arguments are named
## Code from Hadley Wickham's pryr package
all_named = function(x) {
  if (length(x) == 0) return(TRUE)
  !is.null(names(x)) && all(names(x) != "")
}


## Coerce to an environment
## Code from Hadley Wickham's pryr package
to_env = function(x, quiet = FALSE) {
  if (is.environment(x)) {
    x
  } else if (is.list(x)) {
    list2env(x)
  } else if (is.function(x)) {
    environment(x)
  } else if (length(x) == 1 && is.character(x)) {
    if (!quiet) message("Using environment ", x)
    as.environment(x)
  } else if (length(x) == 1 && is.numeric(x) && x > 0) {
    if (!quiet) message("Using environment ", search()[x])
    as.environment(x)
  } else {
    stop("Input can not be coerced to an environment", call. = FALSE)
  }
}


## Create a function from arguments, body code, and environment
## Code from Hadley Wickham's pryr package
make_function = function (args, body, env = parent.frame()) {
  args = as.pairlist(args)
  stopifnot(all_named(args), is.language(body))
  env = to_env(env)
  eval(call("function", args, body), env)
}


## Wrapper to catch warnings, errors, and function call as attributes
## Modified from code written by Martin Morgan and Russell S. Pierce http://stackoverflow.com/a/29465795
catch.fun = function (fun) {
  errorOccurred = FALSE
  fullCall = alist(...=)
  fullCall$fnCall = as.character(substitute(fun))
  make_function(args=fullCall, body=quote({
    fnCall = paste0(fnCall, "(", paste(match.call()[-1L], collapse=", "),")")
    warn = err = NULL
    res = withCallingHandlers(tryCatch(fun(...), error = function(e) {
      err <<- conditionMessage(e)
      errorOccurred <<- TRUE
      NULL
    }), warning = function(w) {
      warn <<- append(warn, conditionMessage(w))
      # invokeRestart("muffleWarning")  # commented to display warnings
    })

    attr(res,"fnCall") = fnCall

    if (errorOccurred) {
      res = "An error occurred in the original function"
    }

    if (is.character(warn)) {
      attr(res,"warningMsg") = warn
    } else {
      attr(res,"warningMsg") = NULL
    }

    if (is.character(err)) {
      attr(res,"errorMsg") = err
    } else {
      attr(res, "errorMsg") = NULL
    }
    return(res)
  }))
}


## Helper functions for catch.fun
## Modified from functions written by Russell S. Pierce
hasWarning = function(x) !is.null(attr(x, "warningMsg"))
hasError = function(x) !is.null(attr(x, "errorMsg"))
isClean = function(x) !(hasError(x) | hasWarning(x))
getWarning = function(x) attr(x, "warningMsg")
getError = function(x) attr(x, "errorMsg")
getFnCall = function(x) attr(x, "fnCall")

