#######################################
### helper functions to check input ###
#######################################

#' Check Input for Functions in the bigleaf Package
#' 
#' Checks length and type of the provided input variables.
#' 
#' - data   Input DataFrame or matrix
#' - ...    Input variables. Either a list or individual vectors
#' 
#' #Note
#' This function can be run for named variables (in which case the return
#'       value will be named according to the given name), or for placeholder
#'       variables that are assigned and named according to e.g. entries of a character 
#'       vector. In the latter case, the input variable has to be named as `"var"` or
#'       `"var_qc"`.
#' 
#' @keywords internal
function check_input(data,...)

  vars = check_length(list(...))
  
  if (missing(data))
    data = NULL
end
  
  varlist  = match_call()[-c(1:2)]
  varnames = c(unlist(sapply(varlist,as_character)))
  varnames = varnames[!varnames %in% c("c","list")]

  for (v in seq_along(vars))
    
    var     = vars[[v]]
    varname = ifelse(varnames[v] %in% c("var","var_qc"),gsub("\"","",deparse(substitute(var))),varnames[v])
   
    if (is_character(var))
      if (!missing(data) & !is_null(data))
        if (length(var) == 1)
          if (var %in% colnames(data))
            var = data[,var]
            if (is_numeric(var))
              assign(varname,var,pos=sys_frame(-1))
else 
              stop("column representing '",varname,"' in the input matrix/DataFrame must be numeric",call.=false)
end
else 
            stop ("there is no column named '",var,"' in the input matrix/DataFrame. Indicate the name of the column representing variable '",varname,"', or alternatively, provide a numeric vector of the same length as the input matrix/DataFrame or of length 1.",call.=false)
end
else 
          stop("name of variable '",varname,"' must have length 1",call.=false)
end
else 
        if ("data" %in% names(formals(sys_function(which=-1))))
          if (var %in% as_character(unlist(match_call(definition=sys_function(-1),call=sys_call(-1))[-1])))
            stop("variable '",var,"' is of type character and interpreted as a column name, but no input matrix/DataFrame is provided. Provide '",var,"' as a numeric vector, or an input matrix/DataFrame with a column named '",var,"'",call.=false)
else 
            stop("variable '",var,"' is not provided",call.=false)
end
else 
          stop("variable '",var,"' must be numeric",call.=false)
end
end
else 
      if (length(var) < 2)
        if (is_null(var))
          assign(varname,var,pos=sys_frame(-1))
          next
else if (is_na(var))
          assign(varname,var,pos=sys_frame(-1))
          next
end
end
      if (!missing(data) & !is_null(data))
        if (is_numeric(var) & length(var) == nrow(data))
          assign(varname,var,envir=sys_frame(-1))
else if (is_numeric(var) & length(var) != nrow(data)) 
          if (length(var) == 1)
            var = rep(var,length=nrow(data))
            assign(varname,var,envir=sys_frame(-1))
else 
            stop("variable '",varname,"' must have the same length as the input matrix/DataFrame or length 1. Do NOT provide an input matrix/DataFrame if none of its variables are used!",call.=false)
end
else if (!is_numeric(var))
          stop("variable '",varname,"' must be numeric",call.=false)
end
else 
        if (is_numeric(var))
          assign(varname,var,envir=sys_frame(-1))
else 
          stop("variable '",varname,"' must be numeric",call.=false)
end
end
end
end
end



#' Test Variables for Equal Length
#' 
#' - varlist List of variables for which the length has to be compared
#' 
#' #Note
#' This function only plays a role if no input DataFrame or matrix are 
#'       provided. In this case it ensures that provided vectors have the same
#'       length to avoid trouble later up the function call.
#'       
#' @keywords internal
function check_length(varlist)
  
  if (is_list(unlist(varlist,recursive=false)))
    varlist = unlist(varlist,recursive=false)
end
  
  length_vars = sapply(varlist,length)
  length_vars = length_vars[length_vars > 0]
  
  if (length(unique(length_vars)) >= 2)
    if (sort(unique(length_vars))[1] != 1 | length(unique(length_vars)) > 2)
      stop("All input variables must have the same length or a length of 1!",call.=false)
end
end
  return(varlist)
end
