##definitions of class unions
##setClassUnion("matrix_Or_numeric_Or_NULL", c("matrix", "numeric", "NULL"))
setClassUnion("matrix_Or_missing", c("matrix", "missing"))
setClassUnion("logical_Or_missing", c("logical", "missing"))
setClassUnion("list_Or_missing", c("list", "missing"))
setClassUnion("character_Or_missing", c("character", "missing"))
setClassUnion("numeric_Or_integer", c("numeric", "integer"))
setClassUnion("numeric_Or_integer_Or_missing", c("numeric", "integer", "missing"))
setClassUnion("integer_Or_numeric_Or_character_Or_missing", c("integer", "numeric", "character", "missing"))