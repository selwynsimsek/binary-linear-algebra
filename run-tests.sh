#!/bin/bash

for lisp in sbcl;
do
$lisp --eval "(time (asdf:test-system 'binary-linear-algebra))" --eval '(uiop:quit)'
done
