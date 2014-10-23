function foo() {

    echo $1
    echo $glvl_var

}

glvl_var=12
foo 45

glvl_var=65
foo 75

