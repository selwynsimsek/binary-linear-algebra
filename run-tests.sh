#!/bin/bash

for lisp in ccl sbcl ecl;
do
$lisp --eval "(time (asdf:test-system 'binary-linear-algebra))" --eval '(uiop:quit)'
done
