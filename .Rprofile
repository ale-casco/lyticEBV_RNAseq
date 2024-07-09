setHook("rstudio.sessionInit", function(newSession) {
  if (newSession)
    rstudioapi::navigateToFile("initiate.R", line = -1L, column = -1L)
}, action = "append")
