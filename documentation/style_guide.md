# Style Guilde

## Shebang

Classic:

    #!/bin/bash

Env-style:

    #!/usr/bin/env bash

* The classic shebang bets everything in one path.
* The env-style shebang looks for the binary in more places and has more chances of running something.

## Comments

    # This is a comment.

* Begin all comments with a capital letter.
* Full stop for sentences.

## Constants

    readonly EXAMPLE_CONSTANT="some constant content"

* Use all caps.
* Use a prefix in order not to duplicate system environment variables.
* Place them at the top of the script.
* Use the ``readonly`` keyword.

## Variables

    example_var="Using underscore"

* Use descriptive name and partial abbreviations if needed.
* Use all lowercase letters and underscores to separate words.
* Use the ``local`` keyword inside functions.

## Commands & exit codes

* Use exit codes to control the script/program execution and to avoid unexpected behavior.


    readonly my_dir="foo"
    ls $my_dir
    ls_return_code=$?
    if [ $ls_code != 0 ] ; then
        echo "List $SOME_DIR directory operation .. Error"
        exit 1
    fi

* 0 stands for 'okay', different than that (1, 2, ..) stands for 'not ok'.

## Functions

Use explicit declaration

    function some_function () {
        return 0
    }
