(uiop:define-package binary-linear-algebra
  (:use #:cl))
(in-package #:binary-linear-algebra)

(defun binary-transpose (matrix)
  "Returns the bt of a binary matrix."
  (destructuring-bind (m n) (array-dimensions matrix)
    (let ((result (bm0 n m)))
      (loop for i from 0 below m
            do (loop for j from 0 below n
                     do (setf (aref result j i) (aref matrix i j))))
      result)))

(define-symbol-macro bt binary-transpose)

(export '(bt binary-transpose))
