#!/bin/bash

for lisp in ccl sbcl;
do
$lisp --eval "(time (asdf:test-system 'binary-linear-algebra))" --eval '(uiop:quit)'
done
