#!/bin/bash

for lisp in ccl;
do
$lisp --eval "(time (asdf:test-system 'binary-linear-algebra))" --eval '(uiop:quit)'
done
