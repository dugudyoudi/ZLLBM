function (AppendFlag  input_command  flag_append)
    set (${input_command} "${${input_command}} ${flag_append}" PARENT_SCOPE)
endfunction (AppendFlag)