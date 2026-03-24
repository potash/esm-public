bind_rows_with_name = function(..., .names_to="name") {
  args = list(...)
  dots = match.call(expand.dots = FALSE)$...
  names = sapply(dots, deparse)
  for(i in 1:length(args)) {
    args[[i]][.names_to] = names[[i]]
  }
  df = bind_rows(args)
  df
}
